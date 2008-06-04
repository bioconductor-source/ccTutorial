\name{barreraExpressionX}
\alias{barreraExpressionX}
\docType{data}
\title{Barrera expression data}
\description{
  Preprocessed expression data for five \code{M. musculus} tissues.
  Samples were hybridized to Affymetrix's mouse430v2 array platform,
  preprocessed using the MAS5 method, as implemented in
  Bioconductor's package \code{affy}.
}
\usage{data(barreraExpressionX)}
\format{
  Preprocessed expression data as an \code{ExpressionSet}
  with 45101 features, 5 samples.
}
\details{
  see the supplement for the code chunk that was used to read in the
  CEL files that are in the \code{expression} directory of the package
  and to preprocess them.
}
%\source{}
%\references{}
\seealso{\code{\link[Biobase]{ExpressionSet-class}},
  \code{\link[affy]{mas5}}}
\examples{
data(barreraExpressionX)
show(barreraExpressionX)
## the raw data is here:
dir(system.file("expression", package="ccTutorial"))
}
\keyword{datasets}
