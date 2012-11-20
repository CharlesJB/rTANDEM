\name{rTResult-class}
\Rdversion{1.1}
\docType{class}
\alias{rTResult-class}

\title{Class \code{"rTResult"}}

\description{
  A rTResult object is designed to contain the information of a X!Tandem
  analysis and allow data mining of this information. 
}

\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("rTResult", ...)}.
}

\section{Slots}{
  \describe{
    \item{\code{result.file}:}{Object of class \code{"character"} ~~ }
    \item{\code{proteins}:}{Object of class \code{"data.table"} ~~ }
    \item{\code{peptides}:}{Object of class \code{"data.table"} ~~ }
    \item{\code{ptm}:}{Object of class \code{"data.table"} ~~ }
    \item{\code{used.parameters}:}{Object of class \code{"data.frame"} ~~ }
    \item{\code{unused.parameters}:}{Object of class \code{"vector"} ~~ }
    \item{\code{sequence.source.paths}:}{Object of class \code{"vector"} ~~ }
    \item{\code{estimated.false.positive}:}{Object of class \code{"integer"} ~~ }
    \item{\code{total.peptides.used}:}{Object of class \code{"integer"} ~~ }
    \item{\code{total.proteins.used}:}{Object of class \code{"integer"} ~~ }
    \item{\code{total.spectra.assigned}:}{Object of class \code{"integer"} ~~ }
    \item{\code{total.spectra.used}:}{Object of class \code{"integer"} ~~ }
    \item{\code{total.unique.assigned}:}{Object of class \code{"integer"} ~~ }
    \item{\code{start.time}:}{Object of class \code{"character"} ~~ }
    \item{\code{xtandem.version}:}{Object of class \code{"character"} ~~ }
    \item{\code{quality.values}:}{Object of class \code{"vector"} ~~ }
    \item{\code{nb.input.models}:}{Object of class \code{"integer"} ~~ }
    \item{\code{nb.input.spectra}:}{Object of class \code{"integer"} ~~ }
    \item{\code{nb.partial.cleavages}:}{Object of class \code{"integer"} ~~ }
    \item{\code{nb.point.mutations}:}{Object of class \code{"integer"} ~~ }
    \item{\code{nb.potential.C.terminii}:}{Object of class \code{"integer"} ~~ }
    \item{\code{nb.potential.N.terminii}:}{Object of class \code{"integer"} ~~ }
    \item{\code{nb.unanticipated.cleavages}:}{Object of class \code{"integer"} ~~ }
    \item{\code{initial.modelling.time.total}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{initial.modelling.time.per.spectrum}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{load.sequence.models.time}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{refinement.per.spectrum.time}:}{Object of class \code{"numeric"} ~~ }
  }
}

\section{Methods}{
No methods defined with class "rTResult" in the signature.
}

\author{
Authors: Frederic Fournier, Charles Joly Beauparlant, Rene Paradis,
Arnaud Droit
Maintainer: Frederic Fournier <frederic.fournier@crchuq.ulaval.ca>,
Charles Joly Beauparlant <charles.joly-beauparlant@crchul.ulaval.ca>
}

\examples{
showClass("rTResult")
}

\keyword{classes}
