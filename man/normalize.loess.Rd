\name{normalize.loess}
\alias{normalize.loess}
\alias{normalize.Cel.loess}
\title{Scale microarray data}
\description{
  normalizes two microarray data to each others.
}
\usage{
normalize.loess(mat, subset=sample(1:(dim(mat)[2]),5000),
                epsilon=10^-2,maxit=1,log.it=T,verbose=T,span=2/3,
                family.loess="symmetric")
normalize.Cel.loess(listcel, ...)
}
\arguments{
  \item{mat}{}
  \item{subset}{}
  \item{epsilon}{}
  \item{maxit}{}
  \item{log.it}{}
  \item{verbose}{}
  \item{span}{}
  \item{family.loess}{}
  \item{listcel}{A \code{Cel.container}}
  \item{...}{any of the options of normalize.loess you would like to modify}
}
\details{
  While a unifired scheme is on its way, please see
  \link{loess.normalize} for a complete descritption...
}
\author{
  Rafael\cr
  (if there is a problem it must be because Laurent butchered his nice code...)
}
\seealso{
  \code{\link{probe.profile}}
}
\keyword{smooth}
