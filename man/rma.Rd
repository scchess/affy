\name{rma}
\alias{rma}
\title{Robust Multi-Array Average expression measure}
\description{
This function converts an \code{AffyBatch} into an \code{exprSet} using
the robust multi-array average (RMA) expression measure.
}
\usage{
rma(object, subset = NULL, verbose = TRUE, normalize = TRUE, bg.correct = TRUE, ...)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{object}{an \code{\link{AffyBatch}}}
  \item{subset}{an indexing vector with the the probe sets do to be
    converted to expression.}
  \item{verbose}{logical value. If \code{TRUE} it writes out some messages.}
  \item{normalize}{logical. If \code{TRUE} it uses quantile normalization.}
  \item{bg.correct}{logical. If \code{TRUE} it uses a model based
    background correction method}
  \item{\dots}{further arguments to be passed to the summary statistic function}
}
\details{
This function creates the expression measure described in Irizarry et al
Biostatistics (2002).
}
\value{
An \code{\link{exprSet}}
}
\references{Irizarry et al. Biostatistics 2002}
\seealso{\code{\link{express}},\code{\link{expresso}} }

\examples{
data(affybatch.example)
eset <- rma(affybatch.example)
}
\keyword{manip}

