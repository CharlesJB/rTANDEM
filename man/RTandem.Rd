\name{RTandem}
\alias{RTandem()}
\alias{RTandem}
\title{
An R encapsulation of X!TANDEM CYCLONE (2010.12.01.1)
}
\description{
The function \code{RTandem()} takes the path of a parameter file as argument and calls X!Tandem on this file.
}
\usage{
RTandem(path)
}

\arguments{
  \item{path}{
    The complete path to the parameter file that X!Tandem will use.
}
}
\value{
  RTandem(path) does not return the results of the computation. Instead
  those results are written to an output file whose location is
  indicated in the parameter file. 
}
\references{
Robertson Craig and Ronald C. Beavis, TANDEM: matching proteins with
mass spectra, Bioinformatics, 2004, 20 1466-7.
http://www.thegpm.org/tandem/ 
}
\author{
Authors: Arnaud Droit and Frederic Fournier
Maintainer: Frederic Fournier <frederic.fournier@crchuq.ulaval.ca>
}

\examples{
# Old X!Tandem call style: we call RTandem(input) on a single
# RTandemParam object.

# We create RTandemParam and from X!Tandem style xml file
# located in the installation folder: 
param <- GetParamFromXML(system.file("extdata/input.xml", package="RTandem"))

# We create a RTandemTaxo object and identify a database for yeast
taxo <- RTandemTaxo()
taxo$peptide.path <- system.file("extdata/taxonomy.xml", package="RTandem")
taxo$peptide.taxon <- "yeast"

# We will adjust those two objects to use one another and to use,
# the path of some data and configuration files located
# in the installation folder:
param$'list path, taxonomy information' <- taxo
param$'list path, default parameters' <- system.file("extdata/default_input.xml", package="RTandem")
param$'spectrum, path' <- system.file("extdata/test_spectra.mgf", package="RTandem")
param$'output, xsl path' <- system.file("extdata/tandem-input-style.xsl", package="RTandem")
param$'output, path' <- paste(getwd(), "output.xml", sep="")

# This is the main command to run RTandem. The statistics of the 
# computation will be written to the console, but the output
# will be written to 'output.xml' in the working directory:
RTandem(input)
}