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
normalize(object, method, ...)

normalize.methods(object)
}
\arguments{
  \item{object}{An \code{AffyBatch}}
  \item{method}{the method used to scale. Different methods will be
    available for the different type of array data.}
  \item{...}{Optional arguments for the normalization funtion. See the
    details below.}
}
\details{
  Generic methods related to array data normalization. Currently defined for
  \code{AffyBatch} objects...

  The method \code{normalize} returns a normalized object, which can be
  used further for analysis or expression values generation.
  An optional parameter \code{method} can be given to specify a method to use.

  A list of \emph{nicknames} for the known methods for an object is returned by
  \code{normalize.methods}. If the nickname of a method is called for
  example "loess", the help page can be accessed by doing \code{normalize?loess}.
  The \emph{nicknames} are stored in
  \code{normalize.AffyBatch.methods}. See the general
  documentation of the package to know more (like how to add your own methods).
}
\author{
  laurent
}
\seealso{
  \code{\link{AffyBatch-class}} and the
  documentations for the respective normalizarion methods.
  }
\examples{
}
\keyword{manip}


