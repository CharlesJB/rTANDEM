\name{tandem}
\alias{tandem}
\alias{rtandem}
\alias{rTANDEM}
\title{
  Calls X!TANDEM CYCLONE (2011.12.01.1) from R objects or xml files
}
\description{
  The function \code{tandem(input)} takes a rTParam object or the path of
  a parameter file as argument and calls X!Tandem on it.
  The function \code{rtandem(data.file, taxon, taxonomy, default.parameters)} is a
  wrapper that can be used to circumvent the need for a rTParam input
  object (or of an xml input file).
}
\usage{
tandem(input)
rtandem(data.file, taxon, taxonomy, default.parameters)
}

\arguments{
  \item{input}{
    A path to a X!Tandem style parameter file or a rTParam object.
  }
  \item{data.file}{
    The path to the file containing the raw data to be analysed (in
    'DTA', 'PKL' or 'MGF' format).
  }
  \item{taxon}{
    A string containing the taxon to be used for the analysis
    (e.g. "yeast" or "Homo sapiens").
  }
  \item{taxonomy}{
    Either a rTTaxo object or the path to a X!Tandem style taxonomy xml
    file.
  }
  \item{default.parameters}{
    Either a rTParam object containing the default parameters to be
    used, or the path to a X!Tandem style default-parameters xml file.
  }
}
\value{
  Both \code{tandem(input)} and \code{rtandem(data.file, taxon,
    taxonomy, default.parameters)} returns the path of the xml output file
  generated.
}
\references{
Robertson Craig and Ronald C. Beavis, TANDEM: matching proteins with
mass spectra, Bioinformatics, 2004, 20 1466-7.
http://www.thegpm.org/tandem/ 
}
\author{
Authors: Frederic Fournier, Charles Joly Beauparlant, Rene Paradis, Arnaud Droit

Maintainer: Frederic Fournier <frederic.fournier@crchuq.ulaval.ca>, Charles Joly Beauparlant <charles.joly-beauparlant@crchul.ulaval.ca>
}

\examples{
# X!Tandem call style: we call tandem(input) on a single
# rTParam object.

# We create rTParam from an X!Tandem xml file
# located in the installation folder: 
param <- GetParamFromXML(system.file("extdata/input.xml", package="rTANDEM"))

# We create a rTTaxo object and identify a database for yeast
taxonomy <- rTTaxo( 
  taxon="yeast",
  format="peptide",
  URL= system.file("extdata/fasta/scd.fasta.pro", package="rTANDEM")
  )

# We will adjust those two objects to use one another and to use
# the path of some data and configuration files located
# in the installation folder:
param <- setParamValue(param, 'list path', 'taxonomy information', taxonomy)
param <- setParamValue(param, 'list path', 'default parameters',
  value=system.file("extdata/default_input.xml", package="rTANDEM"))
param <- setParamValue(param, 'spectrum', 'path',
  value=system.file("extdata/test_spectra.mgf", package="rTANDEM"))
param <- setParamValue(param, 'output', 'xsl path',
  value=system.file("extdata/tandem-input-style.xsl", package="rTANDEM"))
param <- setParamValue(param, 'output', 'path',
  value=paste(getwd(), "output.xml", sep="/"))

# This is the main command to run rTANDEM. The output will be
# written to a file in the working directory and the function
# returns the path to this file.
output.file <- tandem(param)
output.file
}