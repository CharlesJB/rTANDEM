\name{tandem}
\alias{tandem}
\alias{rTANDEM}
\title{
An R encapsulation of X!TANDEM CYCLONE (2011.12.01.1)
}
\description{
The function \code{tandem(input)} takes a rTParam object or the path of a parameter file as argument and calls X!Tandem on it.
}
\usage{
tandem(input)
}

\arguments{
  \item{input}{
    A path to a X!Tandem style parameter file or a rTParam object.
}
}
\value{
  tandem(input) returns the path of the output file generated.
}
\references{
Robertson Craig and Ronald C. Beavis, TANDEM: matching proteins with
mass spectra, Bioinformatics, 2004, 20 1466-7.
http://www.thegpm.org/tandem/ 
}
\author{
Authors: Frederic Fournier, Charles Joly-Beauparlant, Arnaud Droit
Maintainer: Frederic Fournier <frederic.fournier@crchuq.ulaval.ca>
}

\examples{
# Old X!Tandem call style: we call tandem(input) on a single
# rTParam object.

# We create rTParam and from X!Tandem style xml file
# located in the installation folder: 
param <- GetParamFromXML(system.file("extdata/input.xml", package="rTANDEM"))

# We create a rTTaxo object and identify a database for yeast
taxo <- rTTaxo()
taxo$peptide.path <- system.file("extdata/fasta/scd.fasta.pro", package="rTANDEM")
taxo$peptide.taxon <- "yeast"

# We will adjust those two objects to use one another and to use,
# the path of some data and configuration files located
# in the installation folder:
param$'list path, taxonomy information' <- taxo
param$'list path, default parameters' <- system.file("extdata/default_input.xml", package="rTANDEM")
param$'spectrum, path' <- system.file("extdata/test_spectra.mgf", package="rTANDEM")
param$'output, xsl path' <- system.file("extdata/tandem-input-style.xsl", package="rTANDEM")
param$'output, path' <- paste(getwd(), "output.xml", sep="/")

# This is the main command to run rTANDEM. The output will be
# written to a file in the working directory and the function
# returns the path to this file.
output.file <- tandem(param)
}