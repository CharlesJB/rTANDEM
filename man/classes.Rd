\name{classes}
\alias{rTParam}
\alias{rTTaxo}
\alias{rTResult}
\alias{setParamValue}
\title{
  rTANDEM S3 and S4 classes
}
\description{
  The functions \code{rTParam()} and \code{rTTaxo()} instantiate S3
  object of respective class 'rTParam' and 'rTTaxo', while
  \code{setParamValue()} helps set the value of a given parameter in an
  existing rTParam object. 'rTParam' and 'rTTaxo' are
  data.frame specifically structured to be used by the other functions of
  the rTANDEM package. \code{rTResult} is an S4 class structured to hold
  the results of the analysis.
}
\usage{
rTParam()
rTTaxo(taxon=NA, format=NA, URL=NA)
setParamValue(param, category, parameter, value)
}
\arguments{
  \item{param}{
    A rTParam object.
  }
  \item{category}{
    The category of an x!tandem full parameter. For example 'output' in
    "output, histograms" or 'protein' in "protein, cleavage site".
  }
  \item{parameter}{
    The specific parameter of an x!tandem full parameter. For example,
    'histograms' in "output, histograms", or 'cleavage site' in
    "protein, cleavage site".
  }
  \item{value}{
    The value to be assigned to the full parameter of the rTParam
    object.
  }
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