\name{normalize}
\title{Normalize array data}
\docType{methods}
\alias{normalize}
\alias{normalize-methods}
\alias{normalize.methods-methods}
\alias{normalize.methods}
\alias{normalize.methods-methods}
\alias{normalize.AffyBatch.methods}
\description{
  Normalize chips intensity data.
}
\usage{
normalize(object,...)

normalize.methods(object)
}
\arguments{
  \item{object}{An \code{AffyBatch}}
  \item{...}{Optional arguments for the normalization funtion, amond
    which \code{method}. See
    details below.}
}
\details{
  Generic methods related to array data normalization. Currently defined for
  \code{AffyBatch}\ldots

  The method \code{normalize} returns a normalized object, which can be
  used further for analysis or expression values generation.
  An optional parameter \code{method} can be given to specify a method to use.

  A list of \emph{nicknames} for the known methods for an object is returned by
  \code{normalize.methods}. If the nickname of a method is called for
  example "loess", the help page can be accessed by doing \code{normalize?loess}.
  Refer to the vignette describing how to add processing methods for
  further details.
}
\seealso{
  \code{\link{AffyBatch-class}} and the
  documentations for the respective normalization methods.
  }
\examples{
}
\keyword{manip}


