\name{plot.density}
\alias{plot.density}
\alias{plot.density.AffyBatch}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plot.density(x, col = rainbow(dim(x)[2]), ylab = "density", xlab = "x", lwd = rep(2, dim(x)[2]), lty = rep(1, dim(x)[2]), xlim = range(x), ylim = c(0, 1), ...)

plot.density.Affybatch(x, col = rainbow(dim(x)[2]), log=TRUE,
which=c("pm","mm","both"), ylab = "density", xlab = NULL, lwd = rep(2, length(x)), lty = rep(1, length(x)), xlim = NULL, ylim = c(0, 1), ...)

}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{x}{A matrix containing the values to make densities in the columns.}
  \item{object}{A object of clase \code{\link{AffyBatch}}}
  \item{col}{The colors to use fot the different arrays}
  \item{ylab}{a title for the y axis.}
  \item{xlab}{ a title for the x axis.}
  \item{lwd}{A vector of line widths.}
  \item{lty}{A vector of line types.}
  \item{xlim}{the range to be encompassed by the x}
  \item{ylim}{the range to be encompassed by the x}
  \item{\dots}{graphical parameters can be given as arguments to \code{\link{plot}}}
}
\author{Ben Bolstad}
\examples{
data(affybatch.example)
hist(affybatch.example) ##this call plot.density.AffyBatch
}
\keyword{hplot}% at least one, from doc/KEYWORDS
