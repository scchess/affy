\name{normalize.constant}
\alias{normalize.constant}
\title{Scale microarray data}
\description{
  normalize two microarray data to each others.
}
\usage{
normalize.constant(x, refconstant, FUN=mean, na.rm=TRUE)
}
\arguments{
  \item{x}{a vector of intensities on a chip (to normalize to the reference).}
  \item{refconstant}{the constant value used as a reference}
  \item{FUN}{the function applied to x to generate a value}
  \item{na.rm}{an optional parameter for the function FUN}
}
\value{
  a vector of normalized intensities, with an attribute "constant"
  holding the value of the factor used for scaling.
}
\author{
  Laurent\cr
}
\seealso{
  \code{\link{normalize}} to normalize \code{\link{Cel.container}} objects.
}
%\examples{

%
%  ## TO DO...
%}
\keyword{manip}


