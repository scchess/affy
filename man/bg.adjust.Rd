\name{bg.adjust}
\alias{bg.adjust}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Background ajustment (internal function)}
\description{
An internal function to be used as the default with \code{\link{bgc}}.
}
\usage{
bg.adjust(pmmm, n.pts = 2^14)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{pmmm}{a matrix}
  \item{n.pts}{number of points to use in call to \code{density}.}
}
\details{Assumes PMs are a convultion of normal and exponentional. So we
  observe X+Y where X is backround and Y is signal. We use E[Y|X+Y,Y>0]
  as our backround corrected PM.}
\value{a matrix}
\author{Ben Bolstad}
\seealso{\code{\link{bgc}}}

\keyword{manip}% at least one, from doc/KEYWORDS
