\name{generateExprSet-method}
\docType{methods}
\alias{generateExprSet}
\alias{generateExprSet-methods}
\alias{generateExprSet.methods}
\title{generate a set of expression values}
\description{
  Generate a set of expression values from the probe pair
  information. The set of expression is returned as an exprSet object.
}
\usage{
  generateExprSet(x, cdf, method, ...)

  generateExprSet.methods
}
\arguments{
  \item{x}{a \code{\link{Cel.container}} holding the probe level
    informations to generate the expression values.}
  \item{cdf}{a \code{\link{Cdf}} object.}
  \item{method}{the method used to generate the expression value}
  \item{...}{any of the options of normalize.loess you would like to modify}
}
\details{
  The different methods available through this mecanism are stored in
  the variable \code{generateExprSet.methods}.
}
\author{
  Laurent
}
\seealso{
  method \code{generateExprSet} of the class \code{\link{Cel.container}}
}
\examples{
  data(CDF.example)
  data(listcel)

  ids <- c("AFFX-BioC-5_at","AFFX-BioC-3_at")
  
  #eset <- generateExprSet(listcel, CDF.example, "avgdiff", ids=ids)
}
\keyword{manip}



