\name{mm9.g2p}
\alias{mm9.g2p}
\docType{data}
\title{Mapping of Ensembl mm9 genes and reporter match positions}
\description{
  a list holding for each mm9 gene annotated in the Ensembl database, 
  which reporter - if any - have been mapped into its upstream or
  transcribed region.
  Users can create this list themselves and the manuscript contains
  the source code for doing this.
  For convenience sake, this list is also given as a data object.
}
\usage{data(mm9.g2p)}
\format{
  a list with 27434 entries.
}
%\details{  created using biomaRt package in March 2008}
%\source{}
%\references{}
\seealso{\code{\link[Ringo]{features2Probes}}}
\examples{
data(mm9.g2p)
head(mm9.g2p)
## how many reporters are mapped per gene:
if (interactive())
 hist(listLen(mm9.g2p))
}
\keyword{datasets}
