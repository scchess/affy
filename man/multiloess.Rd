\name{multiloess}
\alias{multiloess}
\title{Local Polynomial Regression Fitting}
\description{A modified version of of loess. 
Perform loess for every column of Y, but with the robust weights calculated
using all columns}
\usage{
multiloess(formula, data=NULL, weights, subset, na.action, model = FALSE,
	 span = 0.75, enp.target, degree = 2,
	 normalize = TRUE,
	 family = c("gaussian", "symmetric"),
	 method = c("loess", "model.frame"),
	 control = loess.control(...), ...)
}
\arguments{
See \code{\link[stats]{loess}}.
}
\details{
Please refer to \code{\link[stats]{loess}}.
}
\value{
  See \code{\link[stats]{loess}}.
 }
\references{Astrand, M. (2001) \url{http://www.math.chalmers.se/~magnusaa/maffy/}}
\author{Magnus Astrand}
\seealso{ \code{\link[stats]{loess}}}
\keyword{internal}
