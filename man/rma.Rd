\name{rma}
\alias{rma}
\title{Robust Multi-Array Average expression measure}
\description{
  This function converts an \code{\link{AffyBatch}} into an
  \code{\link[Biobase]{exprSet}} using
the robust multi-array average (RMA) expression measure.
}
\usage{
rma(object,subset=NULL, verbose=TRUE, destructive = FALSE,...)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{object}{an \code{\link{AffyBatch}}}
  \item{subset}{an indexing vector with the the probe sets do to be
    converted to expression. (not currently implemented)}
  \item{verbose}{logical value. If \code{TRUE} it writes out some messages. (Curently ignored)}
  \item{destructive}{logical value. If \code{TRUE} works on the PM matrix in place as much as possible, good for large datasets, but dangerous because it will leave the PM element of the AffyBatch changed after a call to \code{rma}}
  \item{\dots}{further arguments to be passed (not currently implemented - stub for future use)}
}
\details{
This function creates the expression measure described in Irizarry et al
Biostatistics (2003).
}
\value{
An \code{\link{exprSet}}
}

\author{Ben Bolstad \email{bolstad@stat.berkeley.edu}}

\references{Irizarry et al. Biostatistics 2003}
\seealso{\code{\link{express}},\code{\link{expresso}} }

\examples{
data(affybatch.example)
eset <- rma(affybatch.example)
}
\keyword{manip}

