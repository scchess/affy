\name{merge.AffyBatch}
\alias{merge.AffyBatch}
\title{merge two AffyBatch objects}
\description{
  merge two AffyBatch objects into one.
}
\usage{
\method{merge}{AffyBatch}(x, y, annotation = paste(annotation(x),
                 annotation(y)), description = NULL, notes =
                 character(0), ...)
}
\arguments{
  \item{x}{an \code{AffyBatch} object.}
  \item{y}{an \code{AffyBatch} object.}
  \item{annotation}{a \code{character} vector.}
  \item{description}{a \code{characterORmiame}, eventually \code{NULL}.}
  \item{notes}{a \code{character} vector.}
  \item{\dots}{additional arguments.}
}
\details{
  To be done.
}
\value{
  A object if class \code{\link[affy:AffyBatch-class]{AffyBatch}}.
}
\seealso{\code{\link{AffyBatch-class}}}
\keyword{manip}
