\name{conversion}
\alias{GetTaxoFromXML}
\alias{GetParamFromXML}
\alias{GetResultsFromXML}
\alias{WriteParamToXML}
\alias{WriteTaxoToXML}
\title{
Convert X!Tandem xml files to R objects and vice versa  
}
\description{
Functions like \code{GetTaxoFromXML("pathToXML")},
\code{GetParamFromXML("pathToXML")}, \code{GetResultsFromXML(pathToXML)}
creates R object from X!Tandem-style xml files. The functions
\code{WriteParamToXML(paramObject)} and
\code{WriteTaxoToXML(paramObject)} will create an X!Tandem-style xml
file from an R object.
}
\usage{
GetTaxoFromXML(xml.file)
GetParamFromXML(xml.file)
GetResultsFromXML(xml.file)
WriteParamToXML(param, file, embeddedParam=c("write", "skip", "merge"),
embeddedTaxo=c("write","skip") )
WriteTaxoToXML(taxo, file)
}
\arguments{
  \item{xml.file}{
    The path to the xml file that is to be read.
  }
  \item{file}{
    The name of the xml file that is to be created.
  }
  \item{param}{
    An object of class rTParam that will be used to create the
    corresponding xml file.
  }
  \item{taxo}{
    A rTTaxo object whose content will be written to an xml file.
  }
  \item{embeddedParam}{
    The behaviour to adopt if a rTParam object contains an embedded
    rTParam object in the "list path, default parameters" slot. The
    option "merge" will merge the two object together. The option
    "write" will call WriteParamToXML on the embedded rTParam object
    and write it to the the given file name plus suffixe
    "_default_param". It will then replace the embedded object by its
    path in the original object. The option "skip" will just ignore
    this slot.
  }
  \item{embeddedTaxo}{
    The behaviour to adopt if the rTParam object contains an embedded
    rTTaxo object in the "list path, taxonomy information" slot. The
    option "write" will call WriteTaxoToXML on the object and write it
    to the input file name plus suffixe "_taxonomy". It will then
    replace the rTTaxo object by its path in the container rTParam
    object. The option "skip" will just ignore this slot.
  }
}
\value{
  'WriteParamToXML' and 'WriteTaxoToXML' have no return value: they are
  used for their side-effect of creating an xml file.
  'GetTaxoFromXML' returns an object of the S3 class rTTaxo,
  'GetParamFromXML' return an object of the S3 class rTParam, and
  'GetResultsFromXML' returns and object of the S4 class rTResult.
}
\examples{
\dontrun{
# To write a parameter or taxonomy object to a single xml file:
WriteParamToXML(parameter_object, file="parameter_file")
# Produces the file 'parameter_file'.

# To write a parameter object which has an embedded default
# parameter set to two xml files: 
WriteParamToXML(parameter_object, file="param_file",
embeddedParam="write")
# Produces the files 'param_file' and 'param_file_default_param'.

# To write a parameter object that contains an embedded taxonomy to
# two different xml files:
WriteParamToXML(parameter_object, file="parameter_file",
embeddedTaxo="write")
# Produces the files 'parameter_file' and 'parameter_file_taxonomy'.

# To write a taxonomy object to a n xml file:
WriteTaxoToXML(taxonomy_object, file="taxonomy")
# Produces the file 'taxonomy'.

# To read a taxonomy file in R:
taxonomy <- GetTaxoFromXML("taxonomy.xml")
# Read the xml file and create a taxonomy object of class rTTaxo.

# To read a parameter file in R:
param <- GetParamFromXML("parameters.xml")
# Read the xml file and create a parameter object of class rTParam.

# To read a result file in R:
results <- GetResultsFromXML("output.xml")
# Read the output from X!Tandem and creates a R object of class
# rTResult.

}
}