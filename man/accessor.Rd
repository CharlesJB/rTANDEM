\name{accessors}
\alias{GetProteins}
\alias{GetPeptides}
\alias{GetDegeneracy}
\title{
Extract information from rTANDEM result object
}
\description{
The \code{GetProteins}, \code{GetPeptides} and \code{GetDegeneracy}
functions are used to extract information from the rTANDEM result object.
}
\usage{
GetProteins(results, log.expect=0, min.peptides=1L)
GetPeptides(protein.uid, results, expect=1, score=0)
GetDegeneracy(peptide.id, results)
}
\arguments{
  \item{results}{
    An object of the class rTResult that contains the result of an
    rTANDEM or X!Tandem analysis.
  }
  \item{log.expect}{
    X!Tandem provides a score of protein identification that is
    presented in terms of the log of the expect value of the
    identification. This score can be used as a threshold to discard low
    confidence identifications from the protein list.
  }
  \item{expect}{
    The expect value of peptide identification. This statistic can be
    used as a threshold to discard low confidence identifications from
    the peptide list.
  }
  \item{min.peptides}{
    The number of peptides involved in the identification of a given
    protein is computed. This number can be used as a threshold to
    discard identifications based on too few peptides from the protein list.
  }
  \item{protein.uid}{
    The tandem identifier of the protein (a numeric).
  }
  \item{peptide.id}{
    The tandem identifier of the peptide (a character).
  }
  \item{score}{
    The tandem score of the peptide identification. This score can be
    used as a threshold to discard low confidence identifications from
    the peptide list. 
  }
}
\value{
  \code{GetProteins} and  \code{GetDegeneracy} return a data.table
  of proteins. \code{GetPeptides} returns a data.table of peptides with
  their ptm (post-translational modifications). Note that this table is
  generated through a merge of the peptide table and the ptm table:
  hence, if peptides has two ptm, it will occupy to rows in
  the resulting data.table.
}
\examples{
# To show how to use the accessor functions, we need an rTANDEM result.
# We can produce one by running the example from the rTANDEM function,
# and reading it to R with GetResultsFromXML:
# output.file.path <- example(rTANDEM)
\dontshow{
output.file.path <- example(rTANDEM);
}
results <- GetResultsFromXML(output.file.path[[1]])

# To get a data.table of the proteins identified by at least 2 peptides
# and with an expect value of 0.05 or better:
proteins <- GetProteins(results, log.expect=-1.3, min.peptides=2)
proteins[, -c(4,5), with=FALSE] # Colums are removed for better display

# To get a list of the peptides used to identify the first protein
# (YCR012W, uid=576): 
peptides <- GetPeptides(protein.uid="576", results)
peptides

# To get the list of proteins to which proteins a peptide belongs:
# (If a peptide belongs to more than one protein, it should not be
# used for quantification, as a biomarker or a MRM target.)
proteins.of.the.peptide <- GetDegeneracy(peptide.id="169.1.1",
results)
proteins.of.the.peptide[,label]
}