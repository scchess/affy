\name{normalize.quantiles}
\alias{normalize.quantiles}
\title{Quantile Normalization}
\description{Using a normalization based upon quantiles, this function
  normalizes a matrix of probe level intensities.
}
\usage{normalize.quantiles(x)
}
\arguments{
  \item{x}{A matrix of intensities where each column corresponds to a
    chip and each row is a probe.}
}
\details{This method is based upon the concept of a quantile-quantile
  plot extended to n dimensions. No special allowances are made for outliers.
}

\value{a matrix of normalized intensites.}

\author{Ben Bolstad, \email{bolstad@stat.berkeley.edu}}

\seealso{\code{\link{qnrm}}}

\keyword{manip}
