GetProteins <- function(results, log.expect=0L, min.peptides=1L){
  # Get the list of proteins that satisfy certain criteria from a rTResult object.
  # Args:
  #    result: The rTResult object from which the proteins will be extracted.
  #    expect: The maximum log10 of expect value that a protein must have to be listed.
  #            (Expect value of 0.05 correspond to -1.3.)
  #    num.peptides: The minimum number of peptides that must have been identified in a protein
  #            for this protein to be listed.

  ## Dummy declaration to prevent "no visible binding" note when using data.table subset:
  num.peptides=expect.value=NULL
  rm(num.peptides, expect.value)

  prots <-subset(results@proteins, num.peptides >= min.peptides & expect.value < log.expect) 
  setkey(prots, expect.value)
  return(prots)
}

GetPeptides <- function(protein.uid, results, expect=1, score=0){
  # Given a protein uid, get the list of peptides that where used to identify this
  # protein and that satisfy certain criteria.
  # Args:
  #    protein     : The uid of a protein from the result object.
  #    result      : The rTResult object of the experiment.
  #    expect.value: The minimum expectation value for the peptide identification.
  #    tandem.score: The minimum tandem score for the identification.
  # Return:
  #    A data.table of the peptides meeting the criteria and their ptm.

  ## Dummy declaration to prevent "no visible binding" when using data.table subset:
  prot.uid=NULL
  rm(prot.uid)

  setnames(results@ptm, c("type", "at", "modified"), c("ptm.type", "ptm.at", "ptm.modified"))
  peps <- subset(results@peptides, prot.uid==protein.uid)
  return(merge(peps, results@ptm, all.x=TRUE, all.y=FALSE, by="pep.id"))
}

GetDegeneracy <- function(peptide.id, results){
  # Given a peptide.id sequence, get the list of proteins in the result dataset where
  # the sequence of this peptide is found.
  # Args:
  #    peptide: the id of a peptide.
  #    result : The rTResult object of the experiment.
  # Return:
  #    A data.table of the proteins where the sequence of the given peptide is present.

  ## Dummy declaration to prevent "no visible binding" note when using data.table subset:
  pep.id=prot.uid=uid=NULL
  rm(pep.id, prot.uid, uid)
  
  target.seq <- results@peptides[pep.id==peptide.id]$sequence
  prots <- subset(results@peptides, sequence==target.seq, select=prot.uid)
  return(subset(results@proteins, uid %in% prots[[1]]))
}
