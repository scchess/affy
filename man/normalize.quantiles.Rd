\name{normalize.quantiles}
\alias{normalize.quantiles}
\alias{normalize.AffyBatch.quantiles}
\title{Quantile Normalization}
\description{
  Using a normalization based upon quantiles, this function
  normalizes a matrix of probe level intensities.
}
\usage{
  normalize.quantiles(x)
  normalize.AffyBatch.quantiles(abatch, ...)
}
\arguments{
  \item{x}{A matrix of intensities where each column corresponds to a
    chip and each row is a probe.}
  \item{abatch}{An \code{\link{AffyBatch}}}
  \item{...}{extra arguments can be passed (we dont use them yet).} 
}
\details{This method is based upon the concept of a quantile-quantile
  plot extended to n dimensions. No special allowances are made for outliers.
}

\value{
  A normalized \code{AffyBatch}.
}
\references{
  Bolstad, B (2001) \emph{Probe Level Quantile Normalization of High Density
    Oligonucleotide Array Data}. Unpublished manuscript
  \url{http://oz.berkeley.edu/~bolstad/stuff/qnorm.pdf}
  }

\author{Ben Bolstad, \email{bolstad@stat.berkeley.edu}}

\seealso{\code{\link{normalize}}}

\keyword{manip}
