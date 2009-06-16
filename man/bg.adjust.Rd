\name{bg.adjust}
\alias{bg.adjust}
\alias{bg.parameters}
\title{Background adjustment (internal function)}
\description{
An internal function to be used by \code{\link{bg.correct.rma}}.
}
\usage{
bg.adjust(pm, n.pts = 2^14, ...)
bg.parameters(pm, n.pts = 2^14)
}
\arguments{
  \item{pm}{a pm matrix}
  \item{n.pts}{number of points to use in call to \code{density}.}
  \item{\dots}{extra arguments to pass to bg.adjust.}
}
\details{Assumes PMs are a convolution of normal and exponential. So we
  observe X+Y where X is background and Y is signal. \code{bg.adjust}
  returns E[Y|X+Y, Y>0] as our background corrected PM.
  \code{bg.parameters} provides ad hoc estimates of the parameters of the
  normal and exponential distributions.}
\value{a matrix}
\seealso{\code{\link{bg.correct.rma}}}
\keyword{manip}
