\name{maffy.normalise}
\alias{maffy.normalise}
\alias{normalize.contrast}
\title{Normalise Itensities}
\description{Normalises feauture intensities from
  \code{\link{probe.level.object}}s}
\usage{
maffy.normalise(data,subset,verbose=F,span=0.25,family="symmetric",log.it=T)}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{data}{an matrix of intensities.}
  \item{subset}{a vector of indexes describing which probes to use for normalising.}
  \item{verbose}{logical value.}
  \item{span}{See \code{\link[modreg]{loess}}.}
  \item{family}{See \code{\link[modreg]{loess}}.}
  \item{log.it}{logical value.}
}
\details{
Please refer to references.
}

}
\value{
The normalised intensities.
}
\references{
Åstrand, M. (2001) \url{http://www.math.chalmers.se/~magnusaa/maffy/}}
\author{Magnus Åstrand}
\seealso{ \code{\link{maffy.subset}}}
\keyword{internal}





