RTandem <-
  function(input) {
    library(XML)
    ## if input is not a R object of class RTandemParam, we consider it to be a
    ## path to an xml parameter file.
    if(! input %in% ls() || class(input)!=RTandemParam ) {
      input<-.paramFromXML(input)
    }
    RTinput<-.paramFromXML(input)
    RTsexp<-.sexpFromParam(RTinput)
    .Call(RTsexp['param'], RTsexp['peptide'], RTsexp['saps'], RTsexp['mods'], RTsexp['spectrum'])
  }

.sexpFromParam<-
  function(param) {

    ## Get the default parameters
    s_default_param<-param$"list path, default parameters"
    ## if the content of 'list path, default parameters' is not an R object of class RTandemParam,
    ## we consider it to be a path to an xml taxonomy file.
    if(s_default_param %in% ls() && class(eval(as.name(s_default_param)))=="RTandemParam")
      default_param<-eval(as.name(s_default_param))
    else{
      default_param<-.paramFromXML(s_default_param)
    }

    # Merge params and default_param in a single object.
    # In case of conflict, param (i.e. input.xml) will prevails.
    merged_params<-.mergeParams(param,default_param) 
    param_list<- .makeList(merged_params)

    ## if the content of 'list path, taxonomy information' is not an R object of class RTandemTaxo,
    ## we consider it to be a path to an xml parameter file.
    s_taxonomy<-merged_params$"list path, taxonomy information"
    if(s_taxonomy %in% ls()  && class(eval(as.name(s_taxonomy)))=="RTandemTaxo")
      taxonomy<-eval(as.name(s_taxonomy))
    else{
      taxonomy<-.taxoFromXML(s_taxonomy)
    }

    #Create peptide path list (and, saps path list, mod path list, spectrum path list)
    taxa<-strsplit(merged_params$"protein, taxon", split=",[[:space:]]*")
    pep_list<-taxonomy$peptide.path[taxonomy$peptide.taxon %in% taxa]
    sap_list<-taxonomy$saps.path[taxonomy$sap.taxon %in% taxa]
    mod_list<-taxonomy$mods.path[taxonomy$mod.taxon %in% taxa]
    spectrum_list<-taxonomy$spectrum.path[taxonomy$spectrum.taxon %in% taxa]
    
    RTsexp<-list(param=param_list, peptide=pep_list, saps=sap_list, mods=mod_list, spectrum=spectrum_list, check.rows=FALSE)
    return(RTsexp)
  }

.taxoFromXML<-
  function(taxo_file) {
    if (file.access(taxo_file, mode=2) == -1)
      stop(as.character(taxo_file), " cannot be read. Verify that the file exists and that you have the permissions necessary to read it.", call. = TRUE)
    doc<-xmlTreeParse(taxo_file, getDTD=F)   
    root<-xmlRoot(doc)
    taxonomy<-.RTandemTaxo()

#    for(x in xmlChildren
    
  }


.mergeParams<-
  function(param, default_param) {
    # For each elements of default_param, 
    # If param already has this element with a value, we do nothing (values of param are NOT overwritten)
    # Else we add the value to param
    
    for ( col in names(default_param) ) {
      if( col %in% names(param) && !is.na(eval(substitute(param$COL, list(COL=col)))) ){
        next}
      else{
        eval(substitute(param$COL<-default_param$COL, list(COL=col)))}
    }
    return(param)
  }

.makeList<-
  function(merged_params) {
    list_params<-list()
    for (col in names(merged_params) ) {
      val=eval(substitute(merged_params$COL, list(COL=col)))
      if(! is.na(val) ) {
        list_params[[length(list_params)+1]]<-col
        list_params[[length(list_params)+1]]<-val
      }
    }
    return(list_params)
  }                   

.paramFromXML<-
  function(XmlFile) {
    if (file.access(XmlFile, mode=2) == -1)
      stop(as.character(XmlFile), " cannot be read. Verify that the file exists and that you have the permissions necessary to read it.", call. = TRUE)
    doc<-xmlTreeParse(XmlFile, getDTD=F)   
    root<-xmlRoot(doc)
    RTinput<- .RTandemParam()
    
    for (x in xmlChildren(root) ) {
      if(xmlName(x)=="note" && "type" %in% names(xmlAttrs(x)) && xmlAttrs(x)['type']=="input") {

        if( xmlAttrs(x)['label'] %in% names(RTinput) ) {
          if(length(as.character(xmlValue(x))) != 0 ) { # To avoid aving character(0) objects that are a pain to deal with
          eval( substitute(
                           RTinput$label<-value,
                           list(label=as.character(xmlAttrs(x)['label']), value=as.character(xmlValue(x)))
                          ))
        }}
        else {
          eval( substitute(
                           RTinput$label<-value,
                           list(label=as.character(xmlAttrs(x)['label']), value=as.character(xmlValue(x)))
                           ))
          message("\"", as.character(xmlAttrs(x)['label']),
                  "\" is not a parameter described in X!Tandem API. It will still be passed to X!Tandem, but it would be wise to check whether a typo has not slipped in this parameter description (for example, \"spectrum,path\" instead of \"spectrum, path\").")
        }
      }
    } # for x in xmlChildren loop
    return(RTinput)
  } # .inputFromXML2 definition


.RTandemTaxo<-
  function(){
    RTandemTaxo<-data.frame(
                            "peptide"=data.frame("taxon"=NULL,"path"=NULL),
                            "mods"=data.frame("taxon"=NULL, "path"=NULL),
                            "spectrum"=data.frame("taxon"=NULL,"path"=NULL),
                            "saps"=data.frame("taxon"=NULL,"path"=NULL)
                            )
    class(RTandemTaxo)<-"RTandemTaxo"
    return(RTandemTaxo)
  }

.taxoFromXML<-
  function(XmlFile){
    if (file.access(XmlFile, mode=2) == -1)
      stop(as.character(XmlFile), " cannot be read. Verify that the file exists and that you have the permissions necessary to read it.", call. = TRUE)
    
    doc<-xmlTreeParse(XmlFile, getDTD=F)   
    root<-xmlRoot(doc)
    RTtaxo<- .RTandemTaxo()
    for (taxon_node in xmlChildren(root) ) {
      if ( xmlName(taxon_node) == "taxon" ) {
        taxon<-xmlAttrs(taxon_node)['label']
        for( path_node in xmlChildren(taxon_node) ){
          if( xmlName(path_node) =="file" ) {
            format<-xmlAttrs(path_node)["format"]
            path<-xmlAttrs(path_node)["URL"]
            eval(parse(text=sprintf("RTtaxo$%s.path<-c(RTtaxo$%s.path, '%s')", format, format, path)))
            eval(parse(text=sprintf("RTtaxo$%s.taxon<-c(RTtaxo$%s.taxon,'%s')", format, format, taxon)))
          }
        } # end of path loop
      }
    } # end of taxon loop
    return(RTtaxo)
  }# end of .taxoFromXML

.RTandemParam<-
  function() {
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


## for a S4 class of RTandemParam objects
## not ideal because it is difficult to add a new slot on the fly
## setClass("RTandemParam",
##          representation(
##                          "list path, default parameters" = "character",       
##                          "list path, taxonomy information" = "character",     
##                          "output, histogram column width" = "character", 
##                          "output, histograms" = "character",
##                          "output, log path" = "character",
##                          "output, maximum valid expectation value" = "character",
##                          "output, message" = "character",
##                          "output, one sequence copy" = "character",
##                          "output, parameters" = "character",
##                          "output, path" = "character",
##                          "output, path hashing" = "character",
##                          "output, performance" = "character",
##                          "output, proteins" = "character",
##                          "output, results" = "character",
##                          "output, sequence path" = "character",
##                          "output, sort results by" = "character",
##                          "output, sequences" = "character",
##                          "output, spectra" = "character",
##                          "output, xsl path" = "character",
##                          "protein, cleavage C-terminal mass change" = "character",
##                          "protein, cleavage N-terminal mass change" = "character",
##                          "protein, cleavage semi" = "character",
##                          "protein, cleavage site" = "character",
##                          "protein, C-terminal residue modification mass" = "character",
##                          "protein, N-terminal residue modification mass" = "character",
##                          "protein, modified residue mass file" = "character",
##                          "protein, quick acetyl" = "character",
##                          "protein, quick pyrolidone" = "character",
##                          "protein, stP bias" = "character",
##                          "protein, saps" = "character",
##                          "protein, taxon" = "character",
##                          "protein, use annotations" = "character",
##                          "refine, cleavage semi" = "character",
##                          "refine, maximum valid expectation value" = "character",
##                          "refine, modification mass" = "character",
##                          "refine, point mutations" = "character",
##                          "refine, potential modification mass" = "character",
##                          "refine, potential modification motif" = "character",
##                          "refine, potential N-terminus modifications" = "character",
##                          "refine, potential C-terminus modifications" = "character",
##                          "refine, refine" = "character",
##                          "refine" = "character",                             # not in the API, but in the example
##                          "refine, saps" = "character",
##                          "refine, sequence path" = "character",
##                          "refine, spectrum synthesis" = "character",
##                          "refine, tic percent" = "character",
##                          "refine, unanticipated cleavage" = "character",
##                          "refine, use annotations" = "character",
##                          "refine, use potential modifications for full refinement" = "character",
##                          "residue, modification mass" = "character",
##                          "residue, potential modification mass" = "character",
##                          "residue, potential modification motif" = "character",
##                          "scoring, a ions" = "character",
##                          "scoring, b ions" = "character",
##                          "scoring, c ions" = "character",
##                          "scoring, cyclic permutation" = "character",
##                          "scoring, include reverse" = "character",
##                          "scoring, maximum missed cleavage sites" = "character",
##                          "scoring, minimum ion count" = "character",
##                          "scoring, x ions" = "character",
##                          "scoring, y ions" = "character",
##                          "scoring, z ions" = "character",
##                          "spectrum, contrast angle" = "character",
##                          "spectrum, dynamic range" = "character",
##                          "spectrum, fragment mass error" = "character",
##                          "spectrum, fragment mass error units" = "character",
##                          "spectrum, fragment mass type" = "character",
##                          "spectrum, fragment monoisotopic mass error" = "character",
##                          "spectrum, fragment monoisotopic mass error units" = "character",
##                          "spectrum, minimum fragment mz" = "character",
##                          "spectrum, minimum peaks" = "character",
##                          "spectrum, minimum parent m+h" = "character",
##                          "spectrum, neutral loss mass" = "character",
##                          "spectrum, neutral loss window" = "character",
##                          "spectrum, parent monoisotopic mass error minus" = "character",
##                          "spectrum, parent monoisotopic mass error plus" = "character",
##                          "spectrum, parent monoisotopic mass error units" = "character",
##                          "spectrum, parent monoisotopic mass isotope error" = "character",
##                          "spectrum, path" = "character",
##                          "spectrum, path type" = "character",
##                          "spectrum, sequence batch size" = "character",
##                          "spectrum, threads" = "character",
##                          "spectrum, total peaks" = "character",
##                          "spectrum, use neutral loss window" = "character",
##                          "spectrum, use noise suppression" = "character",
##                          "spectrum, use contrast angle" = "character"
##                          ),
##         prototype(
##                         "list path, default parameters" = NULL,       
##                         "list path, taxonomy information" = NULL,     
##                         "output, histogram column width" = NULL, 
##                         "output, histograms" = NULL,
##                         "output, log path" = NULL,
##                         "output, maximum valid expectation value" = NULL,
##                         "output, message" = NULL,
##                         "output, one sequence copy" = NULL,
##                         "output, parameters" = NULL,
##                         "output, path" = NULL,
##                         "output, path hashing" = NULL,
##                         "output, performance" = NULL,
##                         "output, proteins" = NULL,
##                         "output, results" = NULL,
##                         "output, sequence path" = NULL,
##                         "output, sort results by" = NULL,
##                         "output, sequences" = NULL,
##                         "output, spectra" = NULL,
##                         "output, xsl path" = NULL,
##                         "protein, cleavage C-terminal mass change" = NULL,
##                         "protein, cleavage N-terminal mass change" = NULL,
##                         "protein, cleavage semi" = NULL,
##                         "protein, cleavage site" = NULL,
##                         "protein, C-terminal residue modification mass" = NULL,
##                         "protein, N-terminal residue modification mass" = NULL,
##                         "protein, modified residue mass file" = NULL,
##                         "protein, quick acetyl" = NULL,
##                         "protein, quick pyrolidone" = NULL,
##                          "protein, stP bias" = NULL,
##                          "protein, saps" = NULL,
##                          "protein, taxon" = NULL,
##                          "protein, use annotations" = NULL,
##                          "refine, cleavage semi" = NULL,
##                          "refine, maximum valid expectation value" = NULL,
##                          "refine, modification mass" = NULL,
##                          "refine, point mutations" = NULL,
##                          "refine, potential modification mass" = NULL,
##                          "refine, potential modification motif" = NULL,
##                          "refine, potential N-terminus modifications" = NULL,
##                          "refine, potential C-terminus modifications" = NULL,
##                          "refine, refine" = NULL,
##                          "refine" = NULL,
##                          "refine, saps" = NULL,
##                          "refine, sequence path" = NULL,
##                          "refine, spectrum synthesis" = NULL,
##                          "refine, tic percent" = NULL,
##                          "refine, unanticipated cleavage" = NULL,
##                          "refine, use annotations" = NULL,
##                          "refine, use potential modifications for full refinement" = NULL,
##                          "residue, modification mass" = NULL,
##                          "residue, potential modification mass" = NULL,
##                          "residue, potential modification motif" = NULL,
##                          "scoring, a ions" = NULL,
##                          "scoring, b ions" = NULL,
##                          "scoring, c ions" = NULL,
##                          "scoring, cyclic permutation" = NULL,
##                          "scoring, include reverse" = NULL,
##                          "scoring, maximum missed cleavage sites" = NULL,
##                          "scoring, minimum ion count" = NULL,
##                          "scoring, x ions" = NULL,
##                          "scoring, y ions" = NULL,
##                          "scoring, z ions" = NULL,
##                          "spectrum, contrast angle" = NULL,
##                          "spectrum, dynamic range" = NULL,
##                          "spectrum, fragment mass error" = NULL,
##                          "spectrum, fragment mass error units" = NULL,
##                          "spectrum, fragment mass type" = NULL,
##                          "spectrum, fragment monoisotopic mass error" = NULL,
##                          "spectrum, fragment monoisotopic mass error units" = NULL,
##                          "spectrum, minimum fragment mz" = NULL,
##                          "spectrum, minimum peaks" = NULL,
##                          "spectrum, minimum parent m+h" = NULL,
##                          "spectrum, neutral loss mass" = NULL,
##                          "spectrum, neutral loss window" = NULL,
##                          "spectrum, parent monoisotopic mass error minus" = NULL,
##                          "spectrum, parent monoisotopic mass error plus" = NULL,
##                          "spectrum, parent monoisotopic mass error units" = NULL,
##                          "spectrum, parent monoisotopic mass isotope error" = NULL,
##                          "spectrum, path" = NULL,
##                          "spectrum, path type" = NULL,
##                          "spectrum, sequence batch size" = NULL,
##                          "spectrum, threads" = NULL,
##                          "spectrum, total peaks" = NULL,
##                          "spectrum, use neutral loss window" = NULL,
##                          "spectrum, use noise suppression" = NULL,
##                          "spectrum, use contrast angle" = NULL
##                          )
##          )
                         
## .inputFromXML1 <-
##   function(input) {
##     doc<-xmlTreeParse(input, getDTD=F)   
##     root<-xmlRoot(doc)
    
##     RTinput<-xmlApply(root, function(x) {
##       if(xmlName(x)=="note" && "type" %in% names(xmlAttrs(x)) && xmlAttrs(x)['type']=="input")
##         c( label=as.character(xmlAttrs(x)['label']), value=as.character(xmlValue(x)) )
##     })
    
##     RTinput<-unname(RTinput)                    #all elements were named 'note'
##     RTinput<-RTinput[!sapply(RTinput, is.null)]     #remove all NULL from the list
##     return(RTinput)
##   }


## This function uses S4 class
## .paramFromXML<-
##   function(XmlFile) {
##     if (file.access(XmlFile, mode=2) == -1)
##       stop(as.character(XmlFile), " cannot be read. Verify that the file exists and that you have the permissions necessary to read it.", call. = TRUE)
##     doc<-xmlTreeParse(XmlFile, getDTD=F)   
##     root<-xmlRoot(doc)
##     RTinput<- new("RTandemParam")

##     for (x in xmlChildren(root) ) {
##       if(xmlName(x)=="note" && "type" %in% names(xmlAttrs(x)) && xmlAttrs(x)['type']=="input") {

##         if(xmlAttrs(x)['label'] %in% names(getSlots("RTandemParam"))) {
##           eval( substitute(
##                            RTinput@label<-value,
##                            list(label=as.character(xmlAttrs(x)['label']), value=as.character(xmlValue(x)))
##                           ))
##         }
##         else {
##           slot(RTinput,
##                as.character(xmlAttrs(x)['label']),
##                check=FALSE) <- as.character(xmlValue(x))
##           message("\"", as.character(xmlAttrs(x)['label']),
##                   "\" is not a parameter described in X!Tandem API. It will still be passed to X!Tandem, but it would be wise to check whether a typo has not slipped in this parameter description (for example, \"spectrum,path\" instead of \"spectrum, path\").")
##         }          
##       }
##     } # for x in xmlChildren loop
##     return(RTinput)
##   } # .inputFromXML2 definition
