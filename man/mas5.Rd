\name{mas5}
\alias{mas5}
\title{MAS 5.0 expression measure}
\description{
  This function converts an instance of \code{\link{AffyBatch-class}}
  into an instance of \code{\link[Biobase]{exprSet-class}} using
  our implementation of Affymetrix's MAS 5.0 expression measure.
}
\usage{
mas5(object, normalize = TRUE, sc = 500, analysis = "absolute", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{an instance of \code{\link{AffyBatch-class}}}
  \item{normalize}{logical. If \code{TRUE} scale normalization is used
    after we obtain an instance of \code{\link[Biobase]{exprSet-class}}}
  \item{sc}{Value at which all arrays will be scaled to.}
  \item{analysis}{should we do absolute or comparison analysis, although
    "comparison" is still not implemented.}
  \item{\dots}{other arguments to be passed to \code{\link{expresso}}.}
}
\details{
  This function is a wrapper for  \code{\link{expresso}} and \code{\link{affy.scalevalue.exprSet}}.}
\value{
  \code{\link[Biobase]{exprSet-class}} 
}
\seealso{\code{\link{expresso}},\code{\link{affy.scalevalue.exprSet}}}
\examples{
  data(affybatch.example)
  eset <- mas5(affybatch.example)
}
\keyword{manip}% at least one, from doc/KEYWORDS
