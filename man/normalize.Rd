\name{normalize}
\alias{normalize}
\title{Normalize}
\description{
  Normalize probe level data
}
\usage{
normalize(x, method=c("quantile", "loess","contrasts"),
          span=2/3, choose.subset=T, subset.size=5000, verbose=T, maxit=1)
}
\arguments{
  \item{x}{a \code{\link{probe.level.object}}}
  \item{method}{which normalization method to use, "quantile", "loess"
    or "contrasts".}
  \item{span}{\code{\link[modreg]{loess}} span to use with "loess" and "contrast"
    methods}
  \item{choose.subset}{logical value. For the "loess" and "contrast"
    methods, if \code{TRUE} it will do normalization on a subset and predict.}
  \item{subset.size}{size of subset to use in normalization.}
  \item{verbose}{logical value. If \code{TRUE} it writes out some messages.}
  \item{maxit}{maximum number of iterations to use in method "loess".}
}
\details{
  Most of the work for this function is done by either
  \code{\link{normalize.quantiles}}, \code{\link{maffy.normalise}}, or \code{\link{loess.normalize}}}
\value{A \code{\link{probe.level.object}}}
\author{Magus Åstrand, Ben Bolstad, and Rafael A. Irizarry}
\examples{
    data(plob.example)
    affy.mva.pairs(plob.example)
    Data <- normalize(plob.example)
    affy.mva.pairs(Data)

    ##simulated example
    x <- 2^sweep(matrix(rnorm(8000,0,.1),2000,4),1,rnorm(2000,5),"+")
    x[,1] <- x[,1] + sqrt(x[,1])
    x[,3] <- x[,3]+(x[,3])^.33
    ##make fake probe level object
    y <- list(pm=x[1:1000,],mm=x[1001:2000,]+1,nprob=1000,nchips=4,
              chip.names=c("chip 1","chip 2","chip 3","chip 4"))
    affy.mva.pairs(y)
    y <- normalize(y)
    affy.mva.pairs(y)
}
\seealso{\code{\link{maffy.normalise}},  \code{\link{normalize.quantiles}}, \code{\link{loess.normalize}}}
\keyword{smooth}
\keyword{manip}