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
  normalize.AffyBatch.quantiles(abatch, pmonly=FALSE)
}
\arguments{
  \item{x}{A matrix of intensities where each column corresponds to a
    chip and each row is a probe.}
  \item{abatch}{An \code{\link[affy:AffyBatch-class]{AffyBatch}}}
  \item{pmonly}{logical value. If \code{TRUE} only PMs are used.}
}
\details{This method is based upon the concept of a quantile-quantile
  plot extended to n dimensions. No special allowances are made for
  outliers. If you make use of quantile normalization either through
  \code{\link{rma}}, \code{\link{express}} or \code{\link{expresso}}
  please cite Bolstad et al, Bioinformatics (2003).
}

\value{
  A normalized \code{AffyBatch}.
}
\references{
  Bolstad, B (2001) \emph{Probe Level Quantile Normalization of High Density
    Oligonucleotide Array Data}. Unpublished manuscript
  \url{http://oz.berkeley.edu/~bolstad/stuff/qnorm.pdf}

  Bolstad, B. M., Irizarry R. A., Astrand, M, and Speed, T. P. (2003)
  \emph{A Comparison of Normalization Methods for High Density
    Oligonucleotide Array Data Based on Bias and Variance.}
   Bioinformatics 19(2) ,pp 185-193. \url{http://www.stat.berkeley.edu/~bolstad/normalize/normalize.html}
  }

\author{Ben Bolstad, \email{bolstad@stat.berkeley.edu}}

\seealso{\code{\link{normalize}}}

\keyword{manip}
