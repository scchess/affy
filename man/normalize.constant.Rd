\name{normalize.constant}
\alias{normalize.constant}
\alias{normalize.Cel.container.constant}
\alias{normalize.Plob.constant}
\title{Scale Cel objects in a container}
\description{
  Scale chip objects in a container.
}
\usage{
normalize.Cel.container.constant(container, refindex=1, FUN=mean, na.rm=TRUE)
normalize.Plob.constant(plob, refindex, FUN=mean, na.rm=TRUE)
normalize.constant(x, refconstant, FUN=mean, na.rm=TRUE)
}
\arguments{
  \item{container}{a \code{\link[Biobase]{container-class}} (\code{Cel.container} for
    now.. .any kind of array data container in the future...)}
  \item{plob}{ a \code{\link{Plob-class}}}
  \item{x}{a vector of intensities on a chip (to normalize to the reference).}
  \item{refindex}{the index of the array used as a reference.}
  \item{FUN}{A function generating a value from the intensities on an
    array. Typically \code{mean} or \code{median}.}
  \item{na.rm}{Paramater passed to the function FUN.}
}
\value{
  %A \code{\link[Biobase]{container-class}} of normalized objects.
  Respectively a \code{\link{Cel.container}} of normalized \code{Cel} objects, or a
  \code{\link{Plob}}, or a vector of normalized intensities, with an attribute "constant"
  holding the value of the factor used for scaling.
}
\author{
  L. Gautier <laurent@cbs.dtu.dk>
}
\seealso{
  \code{\link{Cel.container}},
  \code{\link{normalize.constant}}
  }
%\examples{
%  
%}
\keyword{manip}
