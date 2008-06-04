\name{mm9.gene2GO}
\alias{mm9.gene2GO}
\docType{data}
\title{GO annotation of Ensembl mm9 genes}
\description{
  a list holding the GO terms annotated for each mm9 gene
  annotated in Ensembl, as retrieved using biomaRt.
  Users can create this list themselves and the supplement contains the
  source code for doing this.
  For convenience sake, this is also given as a data object.
}
\usage{data(mm9.gene2GO)}
\format{
  a list with 10526 entries.
}
\details{
  created using biomaRt package in March 2008
}
%\source{}
%\references{}
\examples{
data(mm9.gene2GO)
head(mm9.gene2GO)
## how many genes are annotated per GO term:
mean(listLen(mm9.gene2GO))
table(listLen(mm9.gene2GO))
}
\keyword{datasets}
