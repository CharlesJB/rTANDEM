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
    The behaviour to addopt if a rTParam object contains and embedded
    rTParam object in the "list path, default parameters" slot. The
    option "merge" will merge the two object together; the option
    "write" will call WriteParamToXML on the embedded rTParam object and
    write it to the the given file name plus suffixe "_default_param"
    and will replace the embedded object by its path in the original
    object, the option "skip" will just ignore this slot.
  }
  \item{embeddedTaxo}{
    The behaviour to adopt if the rTParam object contains and embedded
    rTTaxo object in the "list path, taxonomy information" slot. The
    option "write" will call WriteTaxoToXML on the object and write it
    to the input file name plus suffixe "_taxonomy" and replace the
    rTTaxo object by its path in the container rTParam object; the
    option "skip" will just ignore this slot of the container rTParam
    object. 
  }
}
\value{
  'WriteParamToXML' and 'WriteTaxoToXML' have no return value: they are
  used for their side-effect of creating an xml file.
  'GetTaxoFromXML' returns an object of the S3 class rTTaxo,
  'GetParamFromXML' return an object of the S3 class rTParam, and
  'GetResultsFromXML' returns and object of the S4 class rTResult.
}