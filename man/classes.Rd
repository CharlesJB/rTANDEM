\name{classes}
\alias{rTParam}
\alias{rTTaxo}
\alias{rTResult}
\title{
  rTANDEM S3 and S4 classes
}
\description{
  The functions \code{rTParam()} and \code{rTTaxo()} instantiate S3
  object of respective class 'rTParam' and 'rTTaxo'. Those are
  data.frame specifically structure to be used by the other functions of
  the rTANDEM package. \code{rTResult} is an S4 class structure to hold
  the results of the analysis.
}
\usage{
rTParam()
rTTaxo(taxon=NA, format=NA, URL=NA)
}
\arguments{
  \item{taxon}{
    An optional string or vector containing the name(s) of the taxa
    corresponding to the database files, for example, 'yeast' or 'Homo sapiens'.
  }
  \item{format}{
    An optional string or vector containing the types of the database
    files: 'peptide', 'saps', 'mods' or 'spectrum'.
  }
  \item{URL}{
    An optional string of vector containing the paths of the database
    files.
  }
}