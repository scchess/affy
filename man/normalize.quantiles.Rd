\name{normalize.quantiles}
\alias{normalize.quantiles}
\alias{normalize.Cel.container.quantiles}
\alias{normalize.Plob.quantiles}
\title{Quantile Normalization}
\description{
  Using a normalization based upon quantiles, this function
  normalizes a matrix of probe level intensities.
}
\usage{
  normalize.quantiles(x)
  normalize.Cel.container.quantiles(listcel)
  normalize.Plob.quantiles(object,...)
}
\arguments{
  \item{x}{A matrix of intensities where each column corresponds to a
    chip and each row is a probe.}
  \item{listcel}{A \code{\link{Cel.container}}}
  \item{object}{A \code{\link{Plob}}}
  \item{...}{extra arguments can be passed (we dont use them yet).}

}
\details{This method is based upon the concept of a quantile-quantile
  plot extended to n dimensions. No special allowances are made for outliers.
}

\value{
  a matrix of normalized intensites, or a \code{Cel.container} of
  normalized arrays, or a \code{Plob}.} 

\references{
  Bolstad, B (2001) \emph{Probe Level Quantile Normalization of High Density
    Oligonucleotide Array Data}. Unpublished manuscript
  \url{http://oz.berkeley.edu/~bolstad/stuff/qnorm.pdf}
  }

\author{Ben Bolstad, \email{bolstad@stat.berkeley.edu}}

\seealso{\code{\link{normalize}}}

\keyword{manip}
