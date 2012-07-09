\name{RTandem-package}
\alias{RTandem-package}
\docType{package}
\title{
An R encapsulation of X!TANDEM CYCLONE (2010.12.01.1)
}
\description{
This package encapsulate X!Tandem in R. As now, the package only allows you to call the function RTandem with a path to a valid xml parameter file and the results of the analysis are written to a file 'outside of R'. But in a near future, the package should be extended to include classes for input and output objects so that a complete workflow for protein identification by ms/ms could be available in R. 
}
\details{
\tabular{ll}{
Package: \tab RTandem\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2011-09-08\cr
License: \tab Artistic License 1.0\cr
}
}
\author{
Authors: Arnaud Droit and Frederic Fournier
Maintainer: Frederic Fournier <frederic.fournier@crchuq.ulaval.ca>
}
\references{
Robertson Craig and Ronald C. Beavis, TANDEM: matching proteins with
mass spectra, Bioinformatics, 2004, 20 1466-7.
http://www.thegpm.org/tandem/ 
}
\keyword{ package }
\examples{
# As now, RTandem is used with a single parameter which is the path to a
# parameter file which contain the location of the input and output. Here 
# we rewrite this file to your working directory with the appropriate
# locations of the input files and format files on your system.

# Rewriting the taxonomy file with correct locations of fasta files.
taxonomy<-paste(getwd(), "/example_taxonomy.xml", sep="")
scd<-system.file("extdata/fasta/scd.fasta.pro", package="RTandem")
scd_1<-system.file("extdata/fasta/scd_1.fasta.pro", package="RTandem")
crap<-system.file("extdata/fasta/crap.fasta.pro", package="RTandem")
text<-readLines(system.file("extdata/taxonomy.xml", package="RTandem"))
text<-stringr::str_replace(text, "../fasta/scd.fasta.pro", scd)
text<-stringr::str_replace(text, "../fasta/scd_1.fasta.pro", scd_1)
text<-stringr::str_replace(text, "../fasta/crap.fasta.pro", crap)
writeLines(text, con=taxonomy)

# Rewriting the style sheet file with the correct location of the css file.
style_xsl<-paste(getwd(), "/example_tandem-style.xsl", sep="")
style_css<-system.file("tandem-style.css")
text<-readLines(system.file("tandem-style.xsl", package="RTandem"))
text<-stringr::str_replace(text, "/tandem/tandem-style.css", style_css)
writeLines(text, con=style_xsl)

# Rewriting the default_parameter file with correct location of style sheets:
default_input<-paste(getwd(),"/example_default_input.xml", sep="")
input_style<-system.file("tandem-input-style.xsl")
text<-readLines(system.file("default_input.xml", package="RTandem"))
text<-stringr::str_replace(text, "tandem-input-style.xsl", input_style)
text<-stringr::str_replace(text, "tandem-style.xsl", style_xsl)
writeLines(text, con=default_input)

# Rewriting the input_file with the correct location of the parameter files:
input<-paste(getwd(), "/example_input.xml", sep="")
spectra<-system.file("extdata/test_spectra.mgf", package="RTandem")
output<-paste(getwd(), "/output.xml", sep="")
text<-readLines(system.file("extdata/input.xml", package="RTandem"))
text<-stringr::str_replace(text, "default_input.xml", default_input)
text<-stringr::str_replace(text, "taxonomy.xml", taxonomy)
text<-stringr::str_replace(text, "test_spectra.mgf", spectra)
text<-stringr::str_replace(text, "output.xml", output)
writeLines(text, input)

# This is the main command to run RTandem. The statistics of the 
# computation will be written to the console, but the output
# will be written in the location specified in the input.xml file.
RTandem(input)
}
