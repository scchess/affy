\name{normalize.contrasts}
\alias{normalize.contrasts}
\alias{normalize.AffyBatch.contrasts}
\title{Scale Cel objects in a container}
\description{
  Scale chip objects in a container.
}
\usage{
%normalize.contrast() ## currenlt maffy.normalize
normalize.AffyBatch.contrasts(abatch,span=2/3,choose.subset=TRUE,subset.size=5000,verbose=TRUE,family="symmetric")
}
\arguments{
  \item{abatch}{ an \code{\link{AffyBatch-class}}}
  \item{span}{parameter to be passed to the function \code{\link[modreg]{loess}}.}
  \item{choose.subset}{}
  \item{subset.size}{}
  \item{verbose}{verbosity flag}
  \item{family}{parameter to be passed to the function \code{\link[modreg]{loess}}.}
}
\value{
  An object of the same class as the one passed.
}
\examples{
}

\seealso{
  \code{\link{maffy.normalize}}
  }
%\examples{
%  
%}
\keyword{manip}
