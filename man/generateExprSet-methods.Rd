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
  generateExprSet(x, cdf, method, bg.correct, ...)

  generateExprSet.methods
}
\arguments{
  \item{x}{a \code{\link{Cel.container}} holding the probe level
    informations to generate the expression values.}
  \item{cdf}{a \code{\link{Cdf}} object.}
  \item{method}{the method used to generate the expression value (see
    section 'details').}
  \item{bg.correct}{the method used for background correction (see
    section 'deatils').}
  \item{...}{any of the options of normalize.loess you would like to modify}
}
\details{
  An extra argument \code{ids=} can be passed. It must be a vector of
  affids. The expression values will only be computed and returned for
  these affyids.
  
  The different methods available through this mecanism can be accessed
  by calling the method \code{generateExprSet.methods} with an object of
  call \code{Cel.container} as an argument.

  In the Affymetrix design, \emph{MM} probes were included to measure
  the noise (or background signal). The original algorithm for
  background correction was to subtract the \emph{MM} signal to the
  \emph{PM} signal. The methods currently inclluded in the package are
  "bg.correct.subtractmm", "bg.correct.pmonly" and "bg.correct.adjust"
  (see documentation for the method \code{express}).

  IMPORTANTRUE NOTE: This a 'development version'. Names used for methods,
  including the name of this method
  itself are likely to change in the near future. The aim is to merge
  more completely the \code{Cel.container} and \code{Cdf}, and the
  \code{Plob} object. This method and the fucntion \code{express} will
  be unified.
}
\author{
  Laurent
}
\seealso{
  method \code{generateExprSet} of the class
  \code{\link{Cel.container}}\cr
  \code{\link{express}} 
}
\examples{
  data(CDF.example)
  data(listcel)

  ids <- c("AFFX-BioC-5_at","AFFX-BioC-3_at")
  
  eset <- generateExprSet(listcel, CDF.example, method="avgdiff", bg.correct="bg.correct.subtractmm", ids=ids)
}
\keyword{manip}



