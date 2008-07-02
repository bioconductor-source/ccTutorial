\name{probeAnno}
\alias{probeAnno}
\docType{data}
\title{probeAnno object holding mapping of reporters to genomic positions}
\description{
  This object of class \code{probeAnno} holds the mapping of the
  reporters on the tiling microarrays to positions in the mouse genome
  (assembly: mm9). 
}
\usage{data(probeAnno)}
\format{
  The format is:
Formal class 'probeAnno' [package "Ringo"]
}
\details{
  The package's supplementary vignette shows how this object is
  constructed from the output of the alignment tool Exonerate.
  For sake of speed and for computers with small amounts of RAM, we also
  provide it as a data object here.
}
\seealso{\code{\link[Ringo]{probeAnno-class}}}
\examples{
\dontrun{data("probeAnno"); show(probeAnno);
         chromosomeNames(probeAnno)}
}
\keyword{datasets}
