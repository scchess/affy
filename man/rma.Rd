\name{rma}
\alias{rma}
\title{Robust Multi-Array Average expression measure}
\description{
  This function converts an \code{\link[affy:AffyBatch-class]{AffyBatch}} into an
  \code{\link[Biobase]{exprSet}} using
the robust multi-array average (RMA) expression measure.
}
\usage{
rma(object,subset=NULL, verbose=TRUE, destructive = FALSE,normalize=TRUE,background=TRUE,bgversion=2,...)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{object}{an \code{\link[affy:AffyBatch-class]{AffyBatch}}}
  \item{subset}{an indexing vector with the the probe sets do to be
    converted to expression. (not currently implemented)}
  \item{verbose}{logical value. If \code{TRUE} it writes out some messages. (Curently ignored)}
  \item{destructive}{logical value. If \code{TRUE} works on the PM matrix in place as much as possible, good for large datasets.}
  \item{normalize}{logical value. If \code{TRUE} normalize data using
    quantile normalization}
  \item{background}{logical value. If \code{TRUE} background correct
    using RMA background correction}
  \item{bgversion}{integer value indicating which RMA background to use
1: use background similar to pure R rma background given in affy version 1.0 - 1.0.2
2: use background similar to pure R rma background given in affy version 1.1 and above}
  \item{\dots}{further arguments to be passed (not currently implemented - stub for future use)}
}
\details{
This function computes the RMA (Robust Multichip Average) expression measure described in Irizarry et al
Biostatistics (2003).

Please note that the default background adjustment method was changed during the lead up
to the bioconductor 1.2 release. This means that this function and
\code{link{expresso}} should give results that directly agree.

}
\value{
An \code{\link[affy:exprSet-class]{exprSet}}
}

\author{Ben Bolstad \email{bolstad@stat.berkeley.edu}}

\references{Irizarry et al. Biostatistics 2003}
\seealso{\code{\link{express}},\code{\link{expresso}} }

\examples{
data(affybatch.example)
eset <- rma(affybatch.example)
}
\keyword{manip}

