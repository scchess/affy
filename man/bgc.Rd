\name{bgc}
\alias{bgc}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Background Correction}
\description{
Background corrects the PM values in an object of class
\code{\link{Plob}} using an 
arbitraty function.
}
\usage{
bgc(object, bg = bg.adjust)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{object}{An object of class \code{\link{Plob}}.}
  \item{bg}{a function that determines how background correction is
    done. The default is \code{\link{bg.adjust}}.}
}
\details{
 The function that performs the actual background correction must take
 as an argument a matrix of dimension 2\code{nprobes(object)} x
 \code{nchips(object)}. The intention is to have the PMs in the upper
 half of the matrix and the MMs in the lower half of the matrix and
 returns a \code{nprobes(object)} x \code{nchips(object)} with the
 background corrected PMs. The \code{\link{Plob}} it returns has the
 background corrected PMs in the \code{pm} slot.
}
\value{
  A \code{\link{Plob}}}
\references{Irizarry et al. (2002)}
\author{Rafael A. Irizarry}

\seealso{ \code{\link{bg.adjust}}}

\examples{
     data(Dilution)
     tmp <- bgc(Dilution)
     plot(pm(Dilution)[,1]),pm(tmp)[,1],log="xy",main="PMs before and
          after background correction")
}
\keyword{manip}
