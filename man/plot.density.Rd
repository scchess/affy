\name{plotDensity}
\alias{plotDensity}
\alias{plotDensity.AffyBatch}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Plot Densities}
\description{
Plots the non-parametric density estimates using values contained in the
columns of a matrix.
}
\usage{
plotDensity(mat, col = rainbow(dim(mat)[2]), ylab = "density", xlab = "x", lwd = rep(2, dim(mat)[2]), lty = rep(1, dim(mat)[2]), xlim = range(mat), ylim = c(0, 1), ...)

plotDensity.AffyBatch(x, col = rainbow(dim(x)[2]), log=TRUE,
which=c("pm","mm","both"), ylab = "density", xlab = NULL, lwd = rep(2, length(x)), lty = rep(1, length(x)), xlim = NULL, ylim = c(0, 1), ...)

}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{mat}{A matrix containing the values to make densities in the columns.}
  \item{x}{A object of clase \code{\link{AffyBatch}}}
  \item{log}{logical value. If \code{TRUE} the log of the intensities in
    the \code{AffyBatch} are plotted.}
  \item{which}{should a histogram of the PMs, MMs, or both be made?}
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
hist(affybatch.example) ##this call plotDensity.AffyBatch
}
\keyword{hplot}% at least one, from doc/KEYWORDS






