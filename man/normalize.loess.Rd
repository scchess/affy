\name{normalize.loess}
\alias{normalize.loess}
\alias{normalize.Cel.container.loess}
\alias{normalize.Plob.loess}
\title{Scale microarray data}
\description{
  normalizes two microarray data to each others.
}
\usage{
normalize.loess(mat, subset=sample(1:(dim(mat)[2]),5000),
                epsilon=10^-2,maxit=1,log.it=TRUE,verbose=TRUE,span=2/3,
                family.loess="symmetric")
normalize.Cel.container.loess(listcel, ...)
normalize.Plob.loess(plob, ...)
}
\arguments{
  \item{mat}{a matrix with columns containing the values of the chips
    to normalize.}
  \item{subset}{a subset of the data to fit a loess to.}
  \item{epsilon}{a tolerance value (supposed to be a small value - used
    as a stopping criterium).}
  \item{maxit}{maximum number of iterations.}
  \item{log.it}{logical. If \code{TRUE} it takes the log2 of \code{mat}}
  \item{verbose}{logical. If \code{TRUE} displays current pair of chip being
    worked on.}
  \item{span}{parameter to be passed the function \code{\link[modreg]{loess}}}
  \item{family.loess}{parameter to be passed the function
    \code{\link[modreg]{loess}}. \code{"gaussian"} or \code{"symmetric"}
    are acceptable values for this parameter.}
  \item{listcel}{a \code{\link{Cel.container}} object}
  \item{plob}{a \code{\link{Plob}} object}
  \item{...}{any of the options of normalize.loess you would like to
    modify (described above).}
}
\details{
  While a unifired scheme is on its way, please see
  \link{loess.normalize} for a complete descritption...
}
\author{
  Rafael\cr
  (if there is a problem it must be because I (LG) butchered his nice
  code... )
}
\seealso{
  \code{\link{normalize}}}
\examples{
     data(Dilution)
     x <- pm(Dilution[,1:3])
     mva.pairs(x)
     x <- normalize.loess(x,subset=1:nrow(x))
     mva.pairs(x)
}
\keyword{smooth}


