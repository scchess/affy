\name{generateExprVal.method.playerout}
\alias{generateExprVal.method.playerout}
\alias{playerout.costfunction}
\title{Generate an expression value from the probes informations}
\description{
  Generate an expression from the probes
}
\usage{
  generateExprVal.method.playerout(probes, weights=FALSE)
}
\arguments{
  \item{probes}{a list of \code{probes} slots from \code{PPSet.container}}
  \item{weights}{Should the resulting weights be returned ?}
}
\value{
  A vector of expression values.
}
\details{
A non-parametrical method to weight each perfect match
  probe in the set and to compute a weighted mean of the perfect match
  values. One will notice this method only makes use of the perfect
  matches. (see function \code{playerout.costfunction} for the cost function).
}
\author{
  Laurent <laurent@cbs.dtu.dk>\cr
  (Thanks to E. Lazaridris for the original playerout code and the
  discussions about it)
}
\references{
  \item Lazardris EN, Sinibaldi D, Bloom G, Mane S, Jove R (2001. to be
  published in Mathematical Biosciences in february 2002).
  A simple method to improve probe set estimates from oligonucleotides
  arrays.
}
\keyword{manip}


