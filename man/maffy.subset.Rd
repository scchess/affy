\name{maffy.subset}
\alias{maffy.subset}
\title{Select Subset}
\description{Select a subset of rows with small rank-range over columns.}
\usage{
maffy.subset(data,subset.size=5000,maxit=100,
             subset.delta=max(round(subset.size/100),25),verbose=F)}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{data}{a matrix}
  \item{subset.size}{desired size of subset}
  \item{maxit}{maximum number of iterations}
  \item{subset.delta}{maximum deviation from subset.size}
  \item{verbose}{logical value.}
}
\details{
Please refer to references.
}

}
\value{ 
A list with component
\item{subset}{the indexes for subset.}
}
\references{Åstrand, M. (2001) \url{http://www.math.chalmers.se/~magnusaa/maffy/}}}
\author{Magnus Åstrand}
\seealso{ \code{\link{maffy.normalise}}}
\keyword{internal}
