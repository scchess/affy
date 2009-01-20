\name{generateExprSet-method}
\docType{methods}
\alias{generateExprSet-methods}
\alias{computeExprSet}
\alias{generateExprSet.methods}
\alias{upDate.generateExprSet.methods}

\title{generate a set of expression values}
\description{
  Generate a set of expression values from the probe pair
  information. The set of expression is returned as an \code{\link[Biobase:class.ExpressionSet]{ExpressionSet}} object.
}
\usage{
  computeExprSet(x, pmcorrect.method, summary.method, ...)

  generateExprSet.methods()

  upDate.generateExprSet.methods(x)

}
\arguments{
  \item{x}{a \code{\link[=AffyBatch-class]{AffyBatch}} holding the probe level
    informations to generate the expression values, for computeExprSet,
    and for upDate.generateExprSet.methods it is a character vector..}
  \item{pmcorrect.method}{the method used to correct PM values (see
    section 'details').}
  \item{summary.method}{the method used to generate the expression value (see
    section 'details').}
  \item{...}{any of the options of the normalization you would like to modify}
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
  "bg.correct.subtractmm", "bg.correct.pmonly" and "bg.correct.adjust".

  To alter the available methods for generating ExprSets use
  upDate.generateExprSet.methods.  

}
\seealso{
  method \code{generateExprSet} of the class
  \code{\link[=AffyBatch-class]{AffyBatch}}\cr
  \code{\link{expresso}} 
}
\examples{
if (require(affydata)) {
  data(Dilution)

  ids <- c( "1000_at","1001_at")

  eset <- computeExprSet(Dilution, pmcorrect.method="pmonly",
                         summary.method="avgdiff", ids=ids)
}
}
\keyword{manip}



