\name{normalize.quantiles}
\alias{normalize.quantiles}
\alias{normalize.Cel.quantiles}
\title{Quantile Normalization}
\description{
  Using a normalization based upon quantiles, this function
  normalizes a matrix of probe level intensities.
}
\usage{
  normalize.quantiles(x)
  normalize.Cel.quantiles(listcel)
}
\arguments{
  \item{x}{A matrix of intensities where each column corresponds to a
    chip and each row is a probe.}
  \item{listcel}{A \code{\link{Cel.container}}}
}
\details{This method is based upon the concept of a quantile-quantile
  plot extended to n dimensions. No special allowances are made for outliers.
}

\value{
  a matrix of normalized intensites or a \code{Cel.container} of
  normalized arrays.}

\author{Ben Bolstad, \email{bolstad@stat.berkeley.edu}}

\seealso{\code{\link{qnrm}}}

\keyword{manip}
