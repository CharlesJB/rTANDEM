tandem <- function(input) {
  # Launch X!Tandem on the dataset described in 'input'
  # Args:
  #  input: either a parameter file in xml format appropriate to run X!Tandem
  #         or an R object of the class rTParam
  #
  # Returns:
  #   Creates output files in the directory specified in the parameter object.	    

  
  
  if (class(input) != "rTParam") {
    input <- GetParamFromXML(input)
  }
  
  if ( class(input$'list path, default parameters') == "rTParam") {
    default_param <- input$'list path, default parameters'
  }
  else{
    default_param <- GetParamFromXML(input$'list path, default parameters')
  }

  # MERGE: we merge the input and default param in a single rTParam object.
  # The input parameters have precedence over default_parameters.
  for (col in names(default_param)) {
    if (col %in% names(input) && !is.na(eval(substitute(input$COL, list(COL=col))))) {
       next
    }
    else{
      eval(substitute(input$COL<-default_param$COL, list(COL=col)))
    }
  }

  if ( class(input$'list path, taxonomy information') == "rTTaxo") {
    taxonomy <- input$'list path, taxonomy information'
  }
  else{
    taxonomy <- GetTaxoFromXML(input$'list path, taxonomy information')
  }

  # VECTORISATION: Create a simple vector alterning keys and values
  # to pass to the C++ code. We remove the taxonomy and default_param slots
  # because they will be passed to the C code in another way
  param <- vector(mode="character")
  for (col in names(input)) {
    if( col == 'list path, default parameters' || col == 'list path, taxonomy information' ){
      next
    }
    else{
      val=eval(substitute(input$COL, list(COL=col)))
      if(!is.na(val)) {
        param[[length(param)+1]] <- col
        param[[length(param)+1]] <- val
      }
    }
  }

  taxa <- strsplit(input$"protein, taxon", split=",[[:space:]]*")
  peps <- taxonomy$peptide.path[taxonomy$peptide.taxon %in% taxa]
  saps <- taxonomy$saps.path[taxonomy$sap.taxon %in% taxa]
  mods <- taxonomy$mods.path[taxonomy$mod.taxon %in% taxa]
  spectrum <- taxonomy$spectrum.path[taxonomy$spectrum.taxon %in% taxa]
  
  # the C function tandem takes five character vector containing respectively: 
  #  1) parameter names and parameter values, 2) paths to the peptide databases,
  #  3) paths to the SAPs databases, 4) paths to the mods databases,
  #  5) paths to the spectra databases
  pathName <- .Call("tandem",
        as.vector(param, mode="character"),
        as.vector(peps, mode="character"),
        as.vector(saps, mode="character"),
        as.vector(mods, mode="character"),
        as.vector(spectrum, mode="character"), PACKAGE = "rTANDEM")

  return(pathName)
}

GetTaxoFromXML <- function(xml.file) {
  # Converts an XML file containing an X!Tandem taxonomy to a rTTaxo object
  # Args:
  #   taxo_file: the path to an XML file containing a X!Tandem taxonomy.
  # Returns
  #   an object of class rTTaxo
  
  if (file.access(xml.file, mode=4) == -1) {
    stop(as.character(xml.file), " cannot be read. Verify that the file exists and that you have the permissions necessary to read it.", call. = TRUE)
  }
    
  doc <- xmlTreeParse(xml.file, getDTD=F)   
  root <- xmlRoot(doc)
  RTtaxo <- rTTaxo()
  for (taxon.node in xmlChildren(root)) {
    if (xmlName(taxon.node) == "taxon") {
      taxon<-xmlAttrs(taxon.node)['label']
      for (path.node in xmlChildren(taxon.node)){
        if (xmlName(path.node) == "file") {
          format<-xmlAttrs(path.node)["format"]
          path<-xmlAttrs(path.node)["URL"]
          eval(parse(text=sprintf("RTtaxo$%s.path<-c(RTtaxo$%s.path, '%s')",
                                  format, format, path)))
          eval(parse(text=sprintf("RTtaxo$%s.taxon<-c(RTtaxo$%s.taxon,'%s')",
                                  format, format, taxon)))
        }
      } # end of path loop
    }
  } # end of taxon loop
  return(RTtaxo)
}# end of GetTaxoFromXML

GetParamFromXML<- function(xml.file) {
  # Converts an X!Tandem xml parameter file to a rTParam object
  # Args:
  #   xml.file: the path to an XML file containing X!Tandem parameters
  # Returns:
  #   an object of class rTParam
  
  if (file.access(xml.file, mode=4) == -1)
    stop(as.character(xml.file), " cannot be read. Verify that the file exists and that you have the permissions necessary to read it.", call. = TRUE)
  doc<-xmlTreeParse(xml.file, getDTD=F)   
  root<-xmlRoot(doc)
  RTinput<- rTParam()
    
  for (x in xmlChildren(root) ) {
    if(xmlName(x)=="note" && "type" %in% names(xmlAttrs(x)) && xmlAttrs(x)['type']=="input") {
      if( xmlAttrs(x)['label'] %in% names(RTinput) ) {
        if(length(as.character(xmlValue(x))) != 0 ) { # To avoid aving character(0) objects that are a pain to deal with
          eval(substitute(RTinput$label<-value,
                          list(label=as.character(xmlAttrs(x)['label']),
                               value=as.character(xmlValue(x)))))
        }
      }
      else {
        eval(substitute(RTinput$label<-value,
                        list(label=as.character(xmlAttrs(x)['label']),
                             value=as.character(xmlValue(x)))))
 #       message("\"", as.character(xmlAttrs(x)['label']),
 #               "\" is not a parameter described in X!Tandem API. It will still be passed to X!Tandem, but it would be wise to check whether a typo has not slipped in this parameter description (for example, \"spectrum,path\" instead of \"spectrum, path\").")
      }
    }
  } # for x in xmlChildren loop
  return(RTinput)
} # .inputFromXML2 definition

WriteParamToXML <- function(param, file) {
  # Write an XML file in the X!Tandem param format from a rTParam object
  # Args:
  #    param: a rTParam object containing
  #    file: a string giving the path and name of the output file.
  # Returns:
  #    No return, the function is called for its side effect of creating an xml file.

  if (file.access(file, mode=2) == -1) {
    stop("You can't write to ", as.character(file), ". Check that you have writing permission to this directory.", call. = TRUE)
  }

  
#  saveXML(xmlD
}

rTTaxo <- function() {
  # Constructor method for rTTaxo

  rTTaxo<-data.frame(
                          "peptide"=data.frame("taxon"=NULL,"path"=NULL),
                          "mods"=data.frame("taxon"=NULL, "path"=NULL),
                          "spectrum"=data.frame("taxon"=NULL,"path"=NULL),
                          "saps"=data.frame("taxon"=NULL,"path"=NULL)
                          )
  class(rTTaxo)<-"rTTaxo"
  return(rTTaxo)
}

rTParam <- function() {
  # Constructor method for rTParam
  
  rTParam<-data.frame(
                           "list path, default parameters" = NA,       
                           "list path, taxonomy information" = NA,     
                           "output, histogram column width" = NA, 
                           "output, histograms" = NA,
                           "output, log path" = NA,
                           "output, maximum valid expectation value" = NA,
                           "output, message" = NA,
                           "output, one sequence copy" = NA,
                           "output, parameters" = NA,
                           "output, path" = NA,
                           "output, path hashing" = NA,
                           "output, performance" = NA,
                           "output, proteins" = NA,
                           "output, results" = NA,
                           "output, sequence path" = NA,
                           "output, sort results by" = NA,
                           "output, sequences" = NA,
                           "output, spectra" = NA,
                           "output, xsl path" = NA,
                           "protein, cleavage C-terminal mass change" = NA,
                           "protein, cleavage N-terminal mass change" = NA,
                           "protein, cleavage semi" = NA,
                           "protein, cleavage site" = NA,
                           "protein, C-terminal residue modification mass" = NA,
                           "protein, N-terminal residue modification mass" = NA,
                           "protein, modified residue mass file" = NA,
                           "protein, quick acetyl" = NA,
                           "protein, quick pyrolidone" = NA,
                           "protein, stP bias" = NA,
                           "protein, saps" = NA,
                           "protein, taxon" = NA,
                           "protein, use annotations" = NA,
                           "refine, cleavage semi" = NA,
                           "refine, maximum valid expectation value" = NA,
                           "refine, modification mass" = NA,
                           "refine, point mutations" = NA,
                           "refine, potential modification mass" = NA,
                           "refine, potential modification motif" = NA,
                           "refine, potential N-terminus modifications" = NA,
                           "refine, potential C-terminus modifications" = NA,
                           "refine, refine" = NA,
                           "refine" = NA,
                           "refine, saps" = NA,
                           "refine, sequence path" = NA,
                           "refine, spectrum synthesis" = NA,
                           "refine, tic percent" = NA,
                           "refine, unanticipated cleavage" = NA,
                           "refine, use annotations" = NA,
                           "refine, use potential modifications for full refinement" = NA,
                           "residue, modification mass" = NA,
                           "residue, potential modification mass" = NA,
                           "residue, potential modification motif" = NA,
                           "scoring, a ions" = NA,
                           "scoring, b ions" = NA,
                           "scoring, c ions" = NA,
                           "scoring, cyclic permutation" = NA,
                           "scoring, include reverse" = NA,
                           "scoring, maximum missed cleavage sites" = NA,
                           "scoring, minimum ion count" = NA,
                           "scoring, x ions" = NA,
                           "scoring, y ions" = NA,
                           "scoring, z ions" = NA,
                           "spectrum, contrast angle" = NA,
                           "spectrum, dynamic range" = NA,
                           "spectrum, fragment mass error" = NA,
                           "spectrum, fragment mass error units" = NA,
                           "spectrum, fragment mass type" = NA,
                           "spectrum, fragment monoisotopic mass error" = NA,
                           "spectrum, fragment monoisotopic mass error units" = NA,
                           "spectrum, minimum fragment mz" = NA,
                           "spectrum, minimum peaks" = NA,
                           "spectrum, minimum parent m+h" = NA,
                           "spectrum, neutral loss mass" = NA,
                           "spectrum, neutral loss window" = NA,
                           "spectrum, parent monoisotopic mass error minus" = NA,
                           "spectrum, parent monoisotopic mass error plus" = NA,
                           "spectrum, parent monoisotopic mass error units" = NA,
                           "spectrum, parent monoisotopic mass isotope error" = NA,
                           "spectrum, path" = NA,
                           "spectrum, path type" = NA,
                           "spectrum, sequence batch size" = NA,
                           "spectrum, threads" = NA,
                           "spectrum, total peaks" = NA,
                           "spectrum, use neutral loss window" = NA,
                           "spectrum, use noise suppression" = NA,
                           "spectrum, use contrast angle" = NA,
                           check.names=FALSE
                           )
  class(rTParam)<-"rTParam"
  return(rTParam)
}

