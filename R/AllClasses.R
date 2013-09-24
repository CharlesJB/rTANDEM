setClass("rTResult",
         representation=representation(
           result.file                        ="character",
           proteins                           ="data.table",
           peptides                           ="data.table",
           ptm                                ="data.table",
           used.parameters                    ="data.frame",
           unused.parameters                  ="vector",
           sequence.source.paths              ="vector",
           estimated.false.positive           ="integer",
           total.peptides.used                ="integer",
           total.proteins.used                ="integer",
           total.spectra.assigned             ="integer",
           total.spectra.used                 ="integer",
           total.unique.assigned              ="integer",
           start.time                         ="character",
           xtandem.version                    ="character",
           quality.values                     ="vector",
           nb.input.models                    ="integer",
           nb.input.spectra                   ="integer",
           nb.partial.cleavages               ="integer",
           nb.point.mutations                 ="integer",
           nb.potential.C.terminii            ="integer",
           nb.potential.N.terminii            ="integer",
           nb.unanticipated.cleavages         ="integer",
           initial.modelling.time.total       ="numeric",
           initial.modelling.time.per.spectrum="numeric",
           load.sequence.models.time          ="numeric",
           refinement.per.spectrum.time       ="numeric"
           )
         )


rTTaxo <- function(taxon=NA, format=NA, URL=NA) {
  # Constructor method for rTTaxo
  # Args:
  #    taxon : A taxon for the taxonomy (eg. "yeast") or a vector of taxa.
  #    format: The format of the database (eg. "peptides" or "spectrum") or a vector of formats.
  #    URL   : The path to the database file or a vector of paths.
  # Returns:
  #    A rTTaxo object.
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
### This list all the parameters officially supported in the API
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
                      "spectrum, skyline path" = NA,
                      "spectrum, threads" = NA,
                      "spectrum, total peaks" = NA,
                      "spectrum, use neutral loss window" = NA,
                      "spectrum, use noise suppression" = NA,
                      "spectrum, use contrast angle" = NA,

### This list the parameters that are not officially supported in the API, but are listed as implemented in the TPP documentation.                               
                      "output, http" = NA,
                      "output, sort best scores by" = NA,
                      "output, title" = NA,
                      "protein, cleavage N-terminal limit" = NA,
                      "protein, homolog management" = NA,
                      "protein, use minimal annotations" = NA,
                      "refine, full unanticipated cleavage" = NA,
                      "refine, potential N-terminus modification position limit" = NA,
                      "residue, NG deamidation" = NA,
                      "scoring, algorithm" = NA,
                      "scoring, pluggable scoring" = NA,
                      "spectrum, allowed neutral losses" = NA,
                      "spectrum, check all charges" = NA,
                      "spectrum, homolog error" = NA,
                      "spectrum, maximum parent charge" = NA,
                      "spectrum, use conditioning" = NA,
                      check.names=FALSE
                      )
  class(rTParam)<-c("data.frame", "rTParam")
  return(rTParam)
}

setParamValue <- function(param, category, parameter, value) {
  key <- paste(category, parameter, sep=", ")
  if (is.null(param[[key]])){
    warning("This category/parameter combination is not recognized by rTANDEM. It might be due to a typo, or to the use of an undocumented parameter. After the analysis, check your result's 'used.parameters' and 'unused.parameters' slots to confirm that the parameter was successfully used.")
  }
  param[[key]] <- value
  param
}
      
