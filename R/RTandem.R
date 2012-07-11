RTandem <- function(input) {
  # Launch X!Tandem on the dataset described in 'input'
  # Args:
  #  input: either a parameter file in xml format appropriate to run X!Tandem
  #         or an R object of the class RTandemParam
  #
  # Returns:
  #   Creates output files in the directory specified in the parameter object.	    

  # if input is not a R object of class RTandemParam, we consider it to be a
  # path to an xml parameter file.
  if (!input %in% ls() || class(input) != RTandemParam) {
    input <- GetParamFromXML(input)
  }
  RTsexp <- .sexpFromParam(input)

  # the C function tandem takes five character vector containing respectively: 
  #  1) parameter names and parameter values, 2) paths to the peptide databases,
  #  3) paths to the SAPs databases, 4) paths to the mods databases,
  #  5) paths to the spectra databases
  .Call("tandem",
        as.vector(RTsexp$param, mode="character"),
        as.vector(RTsexp$peptide, mode="character"),
        as.vector(RTsexp$saps, mode="character"),
        as.vector(RTsexp$mods, mode="character"),
        as.vector(RTsexp$spectrum, mode="character"), PACKAGE = "RTandem")
}

GetTaxoFromXML <- function(xml.file) {
  # Converts an XML file containing an X!Tandem taxonomy to a RTandemTaxo object
  # Args:
  #   taxo_file: the path to an XML file containing a X!Tandem taxonomy.
  # Returns
  #   an object of class RTandemTaxo
  
  if (file.access(xml.file, mode=4) == -1) {
    stop(as.character(xml.file), " cannot be read. Verify that the file exists and that you have the permissions necessary to read it.", call. = TRUE)
  }
    
  doc <- xmlTreeParse(xml.file, getDTD=F)   
  root <- xmlRoot(doc)
  RTtaxo <- RTandemTaxo()
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
  # Converts an X!Tandem xml parameter file to a RTandemParam object
  # Args:
  #   xml.file: the path to an XML file containing X!Tandem parameters
  # Returns:
  #   an object of class RTandemParam
  
  if (file.access(xml.file, mode=4) == -1)
    stop(as.character(xml.file), " cannot be read. Verify that the file exists and that you have the permissions necessary to read it.", call. = TRUE)
  doc<-xmlTreeParse(xml.file, getDTD=F)   
  root<-xmlRoot(doc)
  RTinput<- RTandemParam()
    
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
        message("\"", as.character(xmlAttrs(x)['label']),
                "\" is not a parameter described in X!Tandem API. It will still be passed to X!Tandem, but it would be wise to check whether a typo has not slipped in this parameter description (for example, \"spectrum,path\" instead of \"spectrum, path\").")
      }
    }
  } # for x in xmlChildren loop
  return(RTinput)
} # .inputFromXML2 definition

WriteParamToXML <- function(param, file) {
  # Write an XML file in the X!Tandem param format from a RTandemParam object
  # Args:
  #    param: a RTandemParam object containing
  #    file: a string giving the path and name of the output file.
  # Returns:
  #    No return, the function is called for its side effect of creating an xml file.

  if (file.access(file, mode=2) == -1) {
    stop("You can't write to ", as.character(file), ". Check that you have writing permission to this directory.", call. = TRUE)
  }

  
#  saveXML(xmlD
}

.sexpFromParam <- function(param) {
  # Helper function that creates a list of the five character vectors that will be passed
  # to the C function "tandem"
  # Args:
  #   param: the RTandemParam object that contains the parameters to be passed to X!Tandem.
  # Returns:
  #   A named list of the five parameters (param, peptide, saps, mods and spectrum) that we
  #   will pass to the C function "tandem"

  s.default.param <- param$"list path, default parameters"
  # if the content of 'list path, default parameters' is not an R object of class RTandemParam,
  # we consider it to be a path to an xml taxonomy file and we create the R object
  if (s.default.param %in% ls() && class(eval(as.name(s.default.param))) == "RTandemParam") {
    default_param <- eval(as.name(s.default.param))
  }
  else{
    default_param <- GetParamFromXML(s.default.param)
  }
  # Merge params and default_param in a single object.
  # In case of conflict, param (i.e. input.xml) will prevails.
  merged.param <- .mergeParams(param,default_param) 
  param.vec <- .makeVector(merged.param)
  # If the content of 'list path, taxonomy information' is not an R object of class RTandemTaxo,
  # we consider it to be a path to an xml parameter file and we create the R object.
  s.taxonomy <- merged.param$"list path, taxonomy information"
  if (s.taxonomy %in% ls() && class(eval(as.name(s.taxonomy))) == "RTandemTaxo") {
      taxonomy <- eval(as.name(s.taxonomy))
    }
  else{
      taxonomy <- GetTaxoFromXML(s.taxonomy)
    }

  # Create peptide path list (and saps path list, mod path list, spectrum path list, ...)
  taxa <- strsplit(merged.param$"protein, taxon", split=",[[:space:]]*")
  peps <- taxonomy$peptide.path[taxonomy$peptide.taxon %in% taxa]
  saps <- taxonomy$saps.path[taxonomy$sap.taxon %in% taxa]
  mods <- taxonomy$mods.path[taxonomy$mod.taxon %in% taxa]
  spectrum <- taxonomy$spectrum.path[taxonomy$spectrum.taxon %in% taxa]
    
  RTsexp <- list(param=param.vec, peptide=peps, saps=saps, mods=mods, spectrum=spectrum)
  return(RTsexp)
  }

.mergeParams <- function(param, default_param) {
  # Helper function that takes two RTandemParam objects (the input and default_input from
  # X!Tandem) and merges them as one.
  # Args:
  #   param: A RTandemParam object (correspond to X!Tandem input.xml file).
  #   default_param: A RTandemParam object (correspond to X!Tandem default_input.xml file).
  
  for (col in names(default_param)) {
    # If param already has an element, we do nothing: values from param ARE NOT OVERWRITTEN 
    if (col %in% names(param) && !is.na(eval(substitute(param$COL, list(COL=col))))) {
       next
      }
      else{
        eval(substitute(param$COL<-default_param$COL, list(COL=col)))
      }
    }
    return(param)
  }

.makeVector <- function(merged.param) {
  # Helper function that takes a RTandemParam object and return a vector
  # to pass to the C "tandem" function
  # Args:
  #   merged.param: a RTandemParam object
  # Return:
  #   A character vector alternating parameter names and their values.
  
  param.vec <- vector(mode="character")
  for (col in names(merged.param)) {
    val=eval(substitute(merged.param$COL, list(COL=col)))
    if(!is.na(val)) {
      param.vec[[length(param.vec)+1]] <- col
      param.vec[[length(param.vec)+1]] <- val
    }
  }
  return(param.vec)
}                   

RTandemTaxo <- function() {
  # Constructor method for RTandemTaxo

  RTandemTaxo<-data.frame(
                          "peptide"=data.frame("taxon"=NULL,"path"=NULL),
                          "mods"=data.frame("taxon"=NULL, "path"=NULL),
                          "spectrum"=data.frame("taxon"=NULL,"path"=NULL),
                          "saps"=data.frame("taxon"=NULL,"path"=NULL)
                          )
  class(RTandemTaxo)<-"RTandemTaxo"
  return(RTandemTaxo)
}

RTandemParam <- function() {
  # Constructor method for RTandemParam
  
  RTandemParam<-data.frame(
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
  class(RTandemParam)<-"RTandemParam"
  return(RTandemParam)
}

