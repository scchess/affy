\name{normalize.constant}
\alias{normalize.constant}
\title{Scale microarray data}
\description{
  normalizes two microarray data to each others.
}
\usage{
normalize.constant(data, refmean)
}
\arguments{
  \item{data}{a vector of intensities on a chip (to normalize to the reference).}
  \item{refmean}{the mean of the reference array}
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


