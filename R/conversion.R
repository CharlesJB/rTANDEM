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
  saveXML(bioml.node, file=file, prefix='<?xml version="1.0"?>\n')
}

GetResultsFromXML <- function(xml.file) {
  # Parse an X!Tandem result file and create a result object
  # Args:
  #    xml.file: the path to an X!Tandem xml output file.
  # Returns:
  #    an object of the class rTResult

  ## Dummy declaration to prevent "no visible binding" when using data.table subset:
  uid=prot.uid=at=pep.id=NULL
  rm(uid, prot.uid, at, pep.id)
  
if (file.access(xml.file, mode=4) == -1)
  stop(as.character(xml.file), " cannot be read. Verify that the file exists and that you have the permissions necessary to read it.", call. = TRUE)

env1                         <- new.env()

# Result file
env1$results                 <- new('rTResult')
env1$results@used.parameters <- rTParam()
env1$results@result.file     <- xml.file

# Temporary information holders
env1$prot.dt <- data.table("uid"=rep(0L,1000),
                           "expect.value"=rep(0,1000),
                           "label"=rep("",1000),
                           "sequence"=rep("",1000),
                           "file"=rep("",1000),
                           "description"=rep("",1000),
                           "num.peptides"=rep(0L,1000) )

env1$pep.dt  <- data.table("prot.uid"=rep(0L,1000),
                           "pep.id"=rep("",1000),
                           "spectrum.id"=rep(0L,1000),
                           "spectrum.mh"=rep(0,1000),
                           "spectrum.sumI"=rep(0,1000),
                           "spectrum.maxI"=rep(0,1000),
                           "spectrum.fI"=rep(0,1000),
                           "expect.value"=rep(0,1000),
                           "tandem.score"=rep(0,1000),
                           "mh"=rep(0,1000), "delta"=rep(0,1000),
                           "peak.count"=rep(0L,1000),
                           "missed.cleavages"=rep(0L,1000),
                           "start.position"=rep(0L,1000),
                           "end.position"=rep(0L,1000),
                           "sequence"=rep("",1000) )

env1$ptm.dt  <- data.table("pep.id"=rep("",1000),
                           "type"=rep("",1000),
                           "at"=rep(0L,1000),
                           "modified"=rep(0,1000) )

# Booleans and counters for tree navigation
env1$perform.params  <- FALSE
env1$used.params     <- FALSE
env1$unused.params   <- FALSE
env1$prot.key        <- 0L
env1$pep.key         <- 0L
env1$ptm.key         <- 0L
env1$ptm.dt.key      <- 1000L
env1$pep.dt.key      <- 1000L
env1$prot.dt.key     <- 1000L

# node attributes information (to be accessed by the children nodes)
env1$group.attrs     <- NULL
env1$protein.attrs   <- NULL
env1$domain.attrs    <- NULL
env1$content         <- NULL

text <- function(content,...) { ##handler function
  env1$content <- content
}
environment(text) <- env1

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

protein <- function(name, attrs, ...) { #handler function
  env1$protein.attrs <- attrs
  
  env1$prot.key <- env1$prot.key + 1L
  
  if ( env1$prot.dt.key < env1$prot.key) {
    prot2.dt <- data.table("uid"=rep(0L,1000), "expect.value"=rep(0,1000), "label"=rep("",1000), "sequence"=rep("",1000), "file"=rep("",1000), "description"=rep("",1000), "num.peptides"=rep(0L,1000) )
    env1$prot.dt  <- rbindlist(list(env1$prot.dt, prot2.dt))
    env1$prot.dt.key <- env1$prot.dt.key + 1000L
  }
  
  set(env1$prot.dt, env1$prot.key, 1L, as.integer(attrs['uid']))
  set(env1$prot.dt, env1$prot.key, 2L, as.numeric(attrs['expect']))
  set(env1$prot.dt, env1$prot.key, 3L, attrs['label'])
}
environment(protein) <- env1

file <- function(name, attrs, ...) { #handler function
  set(env1$prot.dt, env1$prot.key, 5L,attrs['URL'])
}
environment(file) <- env1

domain <- function(name, attrs, ...) { #handler function
  env1$domain.attrs <- attrs
  env1$pep.key <- env1$pep.key + 1L

  if(! env1$content=="") {
    set(env1$prot.dt, env1$prot.key, 4L, env1$content)
  }

  if( env1$pep.dt.key  < env1$pep.key) {
    pep2.dt  <- data.table("prot.uid"=rep(0L,1000), "pep.id"=rep("",1000), "spectrum.id"=rep(0L,1000), "spectrum.mh"=rep(0.0,1000), "spectrum.sumI"=rep(0.0,1000), "spectrum.maxI"=rep(0.0,1000), "spectrum.fI"=rep(0.0,1000), "expect.value"=rep(0.0,1000), "tandem.score"=rep(0.0,1000), "mh"=rep(0.0,1000), "delta"=rep(0.0,1000), "peak.count"=rep(as.integer(0),1000), "missed.cleavages"=rep(0L,1000), "start.position"=rep(0L,1000), "end.position"=rep(0L,1000), "sequence"=rep("",1000) )
    env1$pep.dt   <- rbindlist(list(env1$pep.dt, pep2.dt))
    env1$pep.dt.key <- env1$pep.dt.key + 1000L
  }
  set(env1$pep.dt, env1$pep.key, 1L, as.integer(env1$protein.attrs['uid']))
  set(env1$pep.dt, env1$pep.key, 2L, env1$domain.attrs['id'] ) 
  set(env1$pep.dt, env1$pep.key, 3L, as.integer(env1$group.attrs['id'])) 
  set(env1$pep.dt, env1$pep.key, 4L, as.numeric(env1$group.attrs['mh'])) 
  set(env1$pep.dt, env1$pep.key, 5L, as.numeric(env1$group.attrs['sumI'])) 
  set(env1$pep.dt, env1$pep.key, 6L, as.numeric(env1$group.attrs['maxI'])) 
  set(env1$pep.dt, env1$pep.key, 7L, as.numeric(env1$group.attrs['fI'])) 
  set(env1$pep.dt, env1$pep.key, 8L, as.numeric(env1$domain.attrs['expect'])) 
  set(env1$pep.dt, env1$pep.key, 9L, as.numeric(env1$domain.attrs['hyperscore'])) 
  set(env1$pep.dt, env1$pep.key, 10L, as.numeric(env1$domain.attrs['mh'])) 
  set(env1$pep.dt, env1$pep.key, 11L, as.numeric(env1$domain.attrs['delta'])) 
  set(env1$pep.dt, env1$pep.key, 12L, as.integer(env1$domain.attrs['peak_count'])) 
  set(env1$pep.dt, env1$pep.key, 13L, as.integer(env1$domain.attrs['missed_cleavages'])) 
  set(env1$pep.dt, env1$pep.key, 14L, as.integer(env1$domain.attrs['start'])) 
  set(env1$pep.dt, env1$pep.key, 15L, as.integer(env1$domain.attrs['end'])) 
  set(env1$pep.dt, env1$pep.key, 16L, as.character(env1$domain.attrs['seq'])) 
}
environment(domain) <- env1

aa <- function(name, attrs, ...) { #handler function
  env1$ptm.key <- env1$ptm.key + 1L
  
  if ( env1$ptm.dt.key < env1$ptm.key ) {
    ptm2.dt  <- data.table("pep.id"=rep("",1000), "type"=rep("",1000), "at"=rep(0L,1000), "modified"=rep(0,1000) )
    env1$ptm.dt   <- rbindlist(list(env1$ptm.dt, ptm2.dt))
    env1$ptm.dt.key <- env1$ptm.dt.key + 1000L
  }
  set(env1$ptm.dt, env1$ptm.key, 1L, env1$domain.attrs['id'])
  set(env1$ptm.dt, env1$ptm.key, 2L, attrs['type'])
  set(env1$ptm.dt, env1$ptm.key, 3L, attrs['at'])
  set(env1$ptm.dt, env1$ptm.key, 4L, attrs['modified'])
}
environment(aa) <- env1

note <- function(x, ...) { #branches function
  label <- as.character(xmlAttrs(x)['label'])
  value <- as.character(xmlValue(x))

  # In protein node
  if (label == "description") {
    # In theory, there could be many different descriptions lines. Here, we only take the last one.
    set(env1$prot.dt, env1$prot.key, 6L, value)
  }
   
  # In group/performance parameters node
  if (env1$perform.params) {
    if( ! is.na(charmatch('list path, sequence source #', label))) {
      label <- 'list path, sequence source'}

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

    # In group/input parameters
  } else if (env1$used.params) {
    attrs <- xmlAttrs(x)
    if ('input' %in% attrs) {
      eval(substitute(env1$results@used.parameters$label<-value,
                      list(label=as.character(attrs['label']),
                           value=as.character(xmlValue(x)))))

    }
    #In group/ unused parameters
  } else if (env1$unused.params) {
    string <- paste(label, value, sep=" = ")
    env1$results@unused.parameters <-
      c(env1$results@unused.parameters, string)
  }
}
environment(note) <- env1


xmlEventParse(xml.file, handlers=list(group=group, protein=protein, file=file, domain=domain, aa=aa, text=text), branches=list(note=note))

# Removing empty rows from prot.dt, pep.dt and ptm.dt
env1$prot.dt <- subset(env1$prot.dt, uid != 0L)
setkey(env1$prot.dt, uid)
env1$prot.dt <- unique(env1$prot.dt)

env1$pep.dt <- subset(env1$pep.dt, prot.uid !=0L)
setkey(env1$pep.dt, prot.uid)

env1$ptm.dt <- subset(env1$ptm.dt, at != 0L)
setkey(env1$ptm.dt, pep.id)

num.peptides.values <- lapply(env1$prot.dt[,uid], function(x){
  nrow(subset(env1$pep.dt, prot.uid==x, select=prot.uid))})
for(i in 1:nrow(env1$prot.dt)){set(env1$prot.dt, i, 7L, num.peptides.values[[i]])}

env1$results@proteins <- env1$prot.dt
env1$results@peptides <- env1$pep.dt
env1$results@ptm      <- env1$ptm.dt

return(env1$results)
}
