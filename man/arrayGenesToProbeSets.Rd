\name{arrayGenesToProbeSets}
\alias{arrayGenesToProbeSets}
\docType{data}
\title{Mapping of Ensembl mm9 genes to Affymetrix Mouse 430v2 probe sets}
\description{
  This is a mapping of Ensembl mm9 genes to Affymetrix Mouse 430v2
  probe sets.
  Onyl Ensembl mm9 genes that mapped at least one probe set on the
  array are included in the list.
  This mapping was created using biomaRt on 27 May 2008.
  Users can create this list themselves and the supplement contains
  the source code for doing this.
  For convenience sake, this object is given as a data object.
}
\usage{data(arrayGenesToProbeSets)}
\format{
  a list with 14643 entries.
}
\details{
  created using biomaRt package in May 2008
}
%\source{}
%\references{}
\examples{
data(arrayGenesToProbeSets)
head(arrayGenesToProbeSets)
table(listLen(arrayGenesToProbeSets))
}
\keyword{datasets}
