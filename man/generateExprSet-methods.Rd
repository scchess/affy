\name{generateExprSet-method}
\docType{methods}
\alias{generateExprSet}
\alias{generateExprSet-methods}
\title{generate a set of expression values}
\description{
  Generate a set of expression values from the probe pair
  information. The set of expression is returned as an exprSet object.
}
\usage{
  generateExprSet(x, cdf, method, ...)
}
\arguments{
  \item{x}{a \code{\link{Cel.container}} holding the probe level
    informations to generate the expression values.}
  \item{cdf}{a \code{\link{Cdf}} object.}
  \item{method}{the method used to generate the expression value}
  \item{...}{any of the options of normalize.loess you would like to modify}
}
\details{
  If the parameter \code{ids} is set to \code{NULL} all the known names
  are taken from \code{cdf@name.levels}.
}
\author{
  Laurent
}
\examples{
  data(CDF.example)
  data(listcel)

  ids <- c("AFFX-BioC-5_at","AFFX-BioC-3_at")
  
  eset <- generateExprSet(listcel, CDF.example, "avgdiff", ids=ids)
}
\keyword{manip}



