\name{normalize.contrasts}
\alias{normalize.contrasts}
\alias{normalize.AffyBatch.contrasts}
\title{Normalize intensities using the contrasts method}
\description{
  Scale chip objects in an \code{\link{AffyBatch-class}}.
}
\usage{
%normalize.contrast() ## currently maffy.normalize
normalize.AffyBatch.contrasts(abatch,span=2/3, choose.subset=TRUE,
                              subset.size=5000, verbose=TRUE,
                              family="symmetric",
                              type=c("together","pmonly","mmonly","separate"))
}
\arguments{
  \item{abatch}{an \code{\link{AffyBatch-class}} object.}
  \item{span}{parameter to be passed to the function
    \code{\link[stats]{loess}}.} 
  \item{choose.subset}{Boolean. Defaults to \code{TRUE}}
  \item{subset.size}{Integer. Number of probesets to use in each subset.}
  \item{verbose}{verbosity flag.}
  \item{family}{parameter to be passed to the function
    \code{\link[stats]{loess}}.} 
  \item{type}{a string specifying how the normalization should be applied.}
}
\value{
  An object of the same class as the one passed.
}
\seealso{
  \code{\link{maffy.normalize}}
}
\keyword{manip}
