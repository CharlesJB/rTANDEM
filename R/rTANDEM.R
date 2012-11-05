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

