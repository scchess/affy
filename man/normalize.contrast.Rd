\name{normalize.contrasts}
\alias{normalize.contrasts}
\alias{normalize.Cel.container.contrasts}
\alias{normalize.Plob.contrasts}
\title{Scale Cel objects in a container}
\description{
  Scale chip objects in a container.
}
\usage{
%normalize.contrast() ## currenlt maffy.normalize
normalize.Cel.container.contrasts(listcel,span=2/3,choose.subset=T,subset.size=5000,verbose=T,family="symmetric")
normalize.Plob.contrasts(object,span=2/3,choose.subset=T,subset.size=5000,verbose=T,family="symmetric")
}
\arguments{
  \item{container}{a \code{\link[Biobase]{container-class}} (\code{Cel.container} for
    now.. .any kind of array data container in the future...)}
  \item{plob}{ a \code{\link{Plob-class}}}
  \item{refindex}{the index of the array used as a reference.}
  \item{FUN}{A function generating a value from the intensities on an
    array. Typically \code{mean} or \code{median}.}
  \item{na.rm}{Paramater passed to the function FUN.}
}
\value{
  An object of the same class as the one passed.
}
\author{
  L. Gautier <laurent@cbs.dtu.dk> (<--- not me... this is a copy/paste
  accident. I think the author is Rafael...)
}
\seealso{
  \code{\link{maffy.normalize}}
  }
%\examples{
%  
%}
\keyword{manip}
