\name{bg.adjust}
\alias{bg.adjust}
\alias{bg.parameters}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Background ajustment (internal function)}
\description{
An internal function to be used as the default with \code{\link{bgc}}.
}
\usage{
bg.adjust(pmmm, n.pts = 2^14)
bg.parameters(pm, mm, n.pts = 2^14)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{pmmm}{a matrix comnbining pm and mm}
  \item{pm}{a pm matrix}
  \item{mm}{a mm matrix}
  \item{n.pts}{number of points to use in call to \code{density}.}
}
\details{Assumes PMs are a convultion of normal and exponentional. So we
  observe X+Y where X is backround and Y is signal. \code{bg.adjust}
    returns E[Y|X+Y,Y>0] as our backround corrected
    PM. \code{bg.parameters} provides adhoc estimates of the parameters
    of the normal and exponential distributions.}
\value{a matrix}
\author{Ben Bolstad}
\seealso{\code{\link{bgc}}}

\keyword{manip}% at least one, from doc/KEYWORDS
