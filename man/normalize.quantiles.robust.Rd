\name{normalize.quantiles.robust}
\alias{normalize.quantiles.robust}
\title{Robusst Quantile Normalization}
\description{Using a normalization based upon quantiles, this function
  normalizes a matrix of probe level intensities. Special allowances are
  made for outliers on the high side.
}
\usage{normalize.quantiles.robust(x,outlier.const=1.5)
}
\arguments{
  \item{x}{A matrix of intensities, columns are chips, rows are probes}
  \item{outlier.const}{Any value greater than outlier.const*IQR + UQ is considered an outlier}
}
\details{This method is based upon
  the concept of a quantile-quantile 
  plot extended to n dimensions. Note that the matrix is of intensities
  not log intensities. The function performs better with raw intensities.
}
\note{This function is still experimental.}
\value{a matrix of normalized intensites}
\author{Ben Bolstad, \email{bolstad@stat.berkeley.edu}}

\seealso{\code{\link{qnrm}}, \code{\link{normalize.quantiles}}}
\keyword{manip}