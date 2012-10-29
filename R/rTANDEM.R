rtandem <- function (data.file, taxon, taxonomy, default.parameters) {
  # Interface to tandem(input). This function createsCreates a basic input file
  #   and launch tandem(input) on this input file.
  # Args:
  #    data.file: The raw data.file that needs to be processed (a mgf or dta file
  #      for example). This correspond to the 'spectrum, path' in the parameter object.
  #    taxon: the species to be searched (e.g. "homo sapiens")
  #    taxonomy: the taxonomy xml or rTTaxo object linking taxon to fasta files.
  #    default.parameters: A X!Tandem style xml parameter file, or an rTParam object
  #      containing all the pertinents parameters
  # Returns:
  #    Creates an output file in the directory specified in the parameters and
  #      returns its path.
  input <- rTParam()
  input$`spectrum, path` <- data.file
  input$`protein, taxon` <- taxon
  input$`list path, taxonomy information` <- taxonomy
  input$`list path, default parameters`   <- default.parameters
  tandem(input)
}

tandem <- function(input) {
  # Launch X!Tandem on the dataset described in 'input'
  # Args:
  #  input: either a parameter file in xml format appropriate to run X!Tandem
  #         or an R object of the class rTParam
  #
  # Returns:
  #   Creates an output file in the directory specified in the parameters and
  #     returns its path.	    
  
  if ( ! "rTParam" %in% class(input) ) {
    input <- GetParamFromXML(input)
  }
  
  if ( 'rTParam' %in% class(input$'list path, default parameters') ) {
    default_param <- input$'list path, default parameters'
  } else{
    default_param <- GetParamFromXML(input$'list path, default parameters')
  }
  
  # MERGE: we merge the input and default param in a single rTParam object.
  # The input parameters have precedence over default_parameters.
  for (col in names(default_param)) {
    if (col %in% names(input) && !is.na(eval(substitute(input$COL, list(COL=col))))) {
      next
    } else{
      eval(substitute(input$COL<-default_param$COL, list(COL=col)))
    }
  }

  if ( 'rTTaxo' %in% class(input$'list path, taxonomy information') ) {
    taxo <- input$'list path, taxonomy information'
   } else{
    taxo <- GetTaxoFromXML(input$'list path, taxonomy information')
  }
  
  # VECTORISATION: Create a simple vector alterning keys and values
  # to pass to the C++ code. We remove the taxonomy and default_param slots
  # because they will be passed to the C code in another way
  param <- vector(mode="character")
  for (col in names(input)) {
    if( col == 'list path, default parameters' || col == 'list path, taxonomy information' ){
      next
    } else{
      val=eval(substitute(input$COL, list(COL=col)))
      if(!is.na(val)) {
        param[[length(param)+1]] <- col
        param[[length(param)+1]] <- val
      }
    }
  }
  
  taxa <- strsplit(input$"protein, taxon", split=",[[:space:]]*")
  peps <- taxo[taxo$taxon %in% taxa & taxo$format=="peptide", 3]
  saps <- taxo[taxo$taxon %in% taxa & taxo$format=="saps", 3]
  mods <- taxo[taxo$taxon %in% taxa & taxo$format=="mods", 3]
  spectrum <- taxo[taxo$taxon %in% taxa & taxo$format=="spectrum", 3]
  
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
  #   xml.file: the path to an XML file containing a X!Tandem taxonomy.
  # Returns
  #   an object of class rTTaxo
  
  if (file.access(xml.file, mode=4) == -1) {
    stop(as.character(xml.file), " cannot be read. Verify that the file exists and that you have the permissions necessary to read it.", call. = TRUE)
  }
  
  taxo.doc <- xmlTreeParse(xml.file, getDTD = FALSE, useInternalNodes = TRUE)   
  length <- length(getNodeSet(taxo.doc, "//file"))
  taxo.root <- xmlRoot(taxo.doc)
  row.count = 1
  taxonomy <- rTTaxo()

  for (taxon.node in xmlChildren(taxo.root)) {
    if (xmlName(taxon.node) == "taxon") {
      taxon<-xmlAttrs(taxon.node)['label']
      for (path.node in xmlChildren(taxon.node)){
        if (xmlName(path.node) == "file") {
          format<-xmlAttrs(path.node)["format"]
          URL<-xmlAttrs(path.node)["URL"]
          taxonomy[row.count,] <- c(taxon, format, URL)
          row.count = row.count + 1
        }
      } # end of path loop
    }
  } # end of taxon loop
  return(taxonomy)
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
      } else{
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

WriteParamToXML <- function(param, file, embeddedParam=c("write","skip","merge"), embeddedTaxo=c("write","skip")) {
  # Write an XML file in the X!Tandem param format from a rTParam object
  # Args:
  #    param: a rTParam object containing
  #    file: a string giving the path and name of the output file.
  #    embeddedParam: if a the input object contains another rTParam object (in the 'list path, default parameter'
  #      slot), the option "merge" will merge the two object together; the option "write" will call WriteParamToXML
  #      on the embedded rTParam object and write it to the the given file name plus suffixe "_default_param" and
  #      will replace the embedded object by its path in the container object, the option "skip" will just ignore
  #      this slot.
  #    embeddedTaxo: if the input object contains a rTTaxo object (in 'list path, taxonomy information'), the
  #      option "write" will call WriteTaxoToXML on the object and write it to the input file name plus suffixe
  #      "_taxonomy" and replace the rTTaxo object by its path in the container rTParam object; the option "skip"
  #      will just ignore this slot of the container rTParam object. 
  # Returns:
  #    No return, the function is called for its side effect of creating an xml file.

#  if (file.access(file, mode=2) == -1) {
#    stop("You can't write to ", as.character(file), ". Check that you have writing permission to this directory.", call. = TRUE)
#  }

  embeddedParam <- match.arg(embeddedParam)
  embeddedTaxo <- match.arg(embeddedTaxo)

  bioml.node <- xmlNode(name="bioml", attrs=NULL, value=NULL)

  for(i in 1:length(param)) {
    
    # Embedded taxonomy file
    if ("rTTaxo" %in% class(param[[i]])){
      if (embeddedTaxo=="skip") {
        next
      }
      if (embeddedTaxo=="write"){
        taxo.file <- paste(file, "_taxonomy", sep="")
        WriteTaxoToXML(param[[i]], taxo.file)
        param[[i]] <- taxo.file
      }
    }
    
    # Embedded parameter file
    if ("rTParam" %in% class(param[[i]])) {
      if (embeddedParam=="skip") {
        next
      }
      if (embeddedParam=="write") {
        param.file <- paste(file, "_default_param", sep="")
        WriteParamToXML(param[[i]], param.file, embeddedParam="skip")
        param[[i]] <- param.file
      }
      if (embeddedParam=="merge") {
        embeddedParam <- param[[i]]
        
        next
      }
    }

    # writing the parameters to file
    if (!is.na(param[[i]])) {
          bioml.node <- addChildren(bioml.node,
                                xmlNode(name="note",
                                        attrs=c(type="input", label=names(param)[[i]]),
                                        value=param[[i]]))
    }
  }
  saveXML(bioml.node, file=file, prefix='<?xml version="1.0"?>\n')
  
}

  

WriteTaxoToXML <- function(taxo, file) {
  # Write an XML file in the X!Tandem taxonomy format from a rTTaxo object.
  # Args:
  #    taxo: a rTTaxo object whose content will be written to an xml file.
  #    file: a string giving the path and name of the file to be created.
  # Returns:
  #    No return, the function is calle for its side effect of creating an xml file.

  bioml.node <- xmlNode(name="bioml", attrs=c(label="x! taxon-to-file matching list"), value=NULL)
  taxa <- levels(as.factor(taxo$taxon))
  
  for (taxon in taxa) {
    taxon.node <- xmlNode(name="taxon", attrs=c(label=taxon))
    for (i in 1:nrow(taxo)){
      print(taxo[[i,1]])
      if( taxo[[i,1]] == taxon) {
        taxon.node <- addChildren(taxon.node,
                                  xmlNode(name="file", attrs=c(format=taxo[[i,2]], URL=taxo[[i,3]])))
      }
    }
    bioml.node <- addChildren(bioml.node, taxon.node)
  }
  print(bioml.node)
  #saveXML(bioml.node, file=file, prefix='<?xml version="1.0"?>\n')
}

rTTaxo <- function(taxon=NA, format=NA, URL=NA) {
  # Constructor method for rTTaxo
  # Args:
  #    taxon : A taxon for the taxonomy (eg. "yeast") or a vector of taxa.
  #    format: The format of the database (eg. "peptides" or "spectrum") or a vector of formats.
  #    URL   : The path to the database file or a vector of paths.
  # Returns:
  #    A rTTaxo object with a pre-assigned data.frame of the given length
  rTTaxo <- data.frame(
                       row.names=NULL,
                       taxon=taxon,
                       format=format,
                       URL=URL
                       )
  class(rTTaxo) <- c('rTTaxo', 'data.frame')
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
  class(rTParam)<-c("data.frame", "rTParam")
  return(rTParam)
}

GetResultsFromXML <- function(xml.file) {
  # Parse an X!Tandem result file and create a result object
  # Args:
  #    xml.file: the path to an X!Tandem xml output file.
  # Returns:
  #    an object of the class rTResult

  if (file.access(xml.file, mode=4) == -1)
    stop(as.character(xml.file), " cannot be read. Verify that the file exists and that you have the permissions necessary to read it.", call. = TRUE)

  env1                     <- new.env()
  env1$results             <- new('rTResult')
  env1$results@used.parameters <- rTParam()
  env1$results@result.file <- xml.file
  env1$group.attrs         <- NULL
  env1$proteins            <- NULL
  env1$perform.params      <- FALSE
  env1$used.params         <- FALSE
  env1$unused.params       <- FALSE
  
  
  group <- function(name, attrs, ...){ ##handler function
    env1$group.attrs <- attrs
    if ("parameters" %in% env1$group.attrs) {
      if ('performance parameters' %in% env1$group.attrs) {
        env1$perform.params <- TRUE
        env1$used.params    <- FALSE
        env1$unused.params  <- FALSE
      } else if ('unused input parameters' %in% env1$group.attrs) {
        env1$perform.params <- FALSE
        env1$used.params    <- FALSE
        env1$unused.params  <- TRUE
      } else if ('input parameters' %in% env1$group.attrs) {
        env1$perform.params <- FALSE
        env1$used.params    <- TRUE
        env1$unused.params  <- FALSE
    } else {
      env1$perform.params <- FALSE
      env1$used.params    <- FALSE
      env1$unused.params  <- FALSE
    }
    }
  }
  environment(group) <- env1

  note <- function(x, ...) { #branche function
    if (env1$perform.params) {
      label <- xmlAttrs(x)['label']
      if( ! is.na(charmatch('list path, sequence source #', label))) {
        label <- 'list path, sequence source'}
      value <- xmlValue(x)
      switch(label,
             'list path, sequence source' = {
               env1$results@sequence.source.paths <-
                 append(env1$results@sequence.source.paths, value)},
             'modelling, estimated false positives' = {
               env1$results@estimated.false.positive <- as.integer(value)},
             'modelling, total peptides used' = {
               env1$results@total.peptides.used <- as.integer(value)},
             'modelling, total proteins used' = {
               env1$results@total.proteins.used <- as.integer(value)},
             'modelling, total spectra assigned' = {
               env1$results@total.spectra.assigned <- as.integer(value)},
             'modelling, total spectra used' = {
               env1$results@total.spectra.used <- as.integer(value)},
             'modelling, total unique assigned' = {
               env1$results@total.unique.assigned <- as.integer(value)},
             'process, start time' = {
               env1$results@start.time <- value},
             'process, version' = {
               env1$results@xtandem.version <- value},
             'quality values' = {
               env1$results@quality.values <- lapply(strsplit(value, "\\s", perl=TRUE), as.integer)},
             'refining, # input models' = {
               env1$results@nb.input.models <- as.integer(value)},
             'refining, # input spectra' = {
               env1$results@nb.input.spectra <- as.integer(value)},
             'refining, # partial cleavage' = {
               env1$results@nb.partial.cleavages <- as.integer(value)},
             'refining, # point mutations' = {
               env1$results@nb.point.mutations <- as.integer(value)},
             'refining, # potential C-terminii' = {
               env1$results@nb.potential.C.terminii <- as.integer(value)},
             'refining, # potential N-terminii' = {
             env1$results@nb.potential.N.terminii <- as.integer(value)},
             'refining, # unanticipated cleavage' = {
               env1$results@nb.unanticipated.cleavages <- as.integer(value)},
             'timing, initial modelling total (sec)' = {
               env1$results@initial.modelling.time.total <- as.numeric(value)},
             'timing, initial modelling/spectrum (sec)' = {
               env1$results@initial.modelling.time.per.spectrum <- as.numeric(value)},
             'timing, load sequence models (sec)' = {
               env1$results@load.sequence.models.time <- as.numeric(value)},
             'timing, refinement/spectrum (sec)' = {
               env1$results@refinement.per.spectrum.time <- as.numeric(value)}
             )
    } else if (env1$used.params) {
      attrs <- xmlAttrs(x)
      if ('input' %in% attrs) {
        eval(substitute(env1$results@used.parameters$label<-value,
                        list(label=as.character(attrs['label']),
                             value=as.character(xmlValue(x)))))
        
      }
    } else if (env1$unused.params) {
      label <- as.character(xmlAttrs(x)['label'])
      value <- as.character(xmlValue(x))
      string <- paste(label, value, sep=" = ")
      env1$results@unused.parameters <-
        append(env1$results@unused.parameters, string)
    }
  }
  environment(note) <- env1
  
  protein <- function(x,...){ #branche function
    spectrum <- c(
                  "id"   = as.integer(env1$group.attrs['id']),
                  "mh"   = as.numeric(env1$group.attrs['mh']),
                  "sumI" = as.numeric(env1$group.attrs['sumI']),
                  "maxI" = as.numeric(env1$group.attrs['maxI']),
                  "fI"   = as.numeric(env1$group.attrs['fI'])
                  )
    protein.attrs <- xmlAttrs(x)
    key <- as.character(protein.attrs['id'])
    pep.node <- xmlElementsByTagName(x, name="peptide")[[1]]
    protein <- list(
                    "uid"          = as.numeric(protein.attrs['uid']),
                    "expect.value" = as.numeric(protein.attrs['expect']),
                    "label"        = protein.attrs['label'],
                    "file"         = xmlAttrs(xmlElementsByTagName(x, name="file")[[1]])['URL'],
                    "description"  = xmlValue(xmlElementsByTagName(x, name="note")[[1]]),
                    "sequence"     = xmlValue(pep.node),
                    "peptides"     = lapply(xmlElementsByTagName(pep.node,name="domain"), function(domain.node) {
#                  "peptides" = xpathApply(pep.node, "./domain", function(domain.node) {
                      domain.attrs <- xmlAttrs(domain.node)
                      peptide <- list(
                         "expect.value" = as.numeric(domain.attrs['expect']),
                         "tandem.score" = as.numeric(domain.attrs['hyperscore']),
                         "mh" = as.numeric(domain.attrs['mh']),
                         "delta" = as.numeric(domain.attrs['delta']),
                         "peak.count" = as.integer(domain.attrs['peak_count']),
                         "missed.cleavage" = as.integer(domain.attrs['missed_cleavages']),
                         "start.position" = as.integer(domain.attrs['start']),
                         "end.position" = as.integer(domain.attrs['end']),
                         "sequence" =  domain.attrs['seq'],
                         "spectrum" = spectrum,
                         "ptm"      = lapply(xmlElementsByTagName(domain.node, name="aa"), function(aa.node){
#                         "ptm" = xpathApply(domain.node, "./aa", function(aa.node){
#                          (xpathApply is 3 sec. slower than lapply(xmlElementsByTagName)
                        #  but it preserves the names and structure of the 'ptm' elements
                           ptm.attrs <- xmlAttrs(aa.node)
                           ptm <- list(
                                       "type" = ptm.attrs['type'],
                                       "position" = as.integer(ptm.attrs['at']),
                                       "mass.change" = as.numeric(ptm.attrs['modified'])
                                       )
                           return(ptm)
                         })
                                      ) # end of peptide <- list
                      return(peptide)
                    })
                    ) # end of protein <- list
    env1$proteins[[key]] <- protein
  }
  environment(protein) <- env1
  
  xmlEventParse(xml.file, handlers=list(group=group), branches=list(protein=protein, note=note))
  
  
  protein.list <- as.list(env1$proteins)
  uid.vector  <- lapply(protein.list, function(x) x$uid)
  duplicated  <- duplicated(uid.vector)
  first.match <- match(uid.vector, uid.vector)
  
  for (i in 1:length(uid.vector)) {
    if (duplicated[i]) {
      protein.list[[first.match[i]]]$peptides <-
        append(protein.list[[first.match[i]]]$peptides, protein.list[[i]]$peptides)
    }
  }
  env1$ results@proteins <- protein.list[!duplicated]
  
  return(env1$results)
}
