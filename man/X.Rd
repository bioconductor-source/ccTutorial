\name{X}
\alias{X}
\docType{data}
\title{ExpressionSet containing ChIP-chip data}
\description{
  Preprocessed ChIP-chexpression data from \emph{M. musculus} brain and
  heart cells.
  H3K4me3 ChIP and input samples were hybridized to a set of four custom
  NimbleGen microarrays with non-overlapping sets of 390k reporters
  each. The results were preprocessed per array type and the result fold
  changes were combinded via \code{rbind} afterwards, resulting in this
  ExpressionSet with 1.5m reporters.
}
\usage{data(X)}
\format{
  Preprocessed ChIP-chip data as an \code{ExpressionSet}
  with 1495582 features, 2 samples.
}
\details{
  This ExpressionSet is created in the vignette, but for purposes of
  speed and for computers with small RAM, we provide it as a data object
  here as well.
}
%\source{}
\seealso{\code{\link[Biobase]{ExpressionSet-class}},
  \code{\link[Ringo]{preprocess}}}
\examples{
data("X"); show(X); pData(X)
}
\keyword{datasets}
