\name{normalize.invariantset}
\alias{normalize.invariantset}
\alias{normalize.Cel.container.invariantset}
\alias{normalize.Plob.invariantset}
\title{Scale Cel objects in a container}
\description{
  Scale chip objects in a container.
}
\usage{
normalize.Cel.container.invariantset(container, f.cdf, prd.td=c(0.003,
0.007), progress=FALSE)

normalize.Plob.invariantset(plob, prd.td=c(0.003, 0.007),
progress=FALSE)

normalize.invariantset(data, ref, prd.td=c(0.003,0.007))
}
\arguments{
  \item{container}{a \code{\link[Biobase]{container-class}} (\code{Cel.container} for
    now.. .any kind of array data container in the future...)}
  \item{plob}{a \code{Plob} (surprise)}
  \item{data}{a vector of intensities on a chip (to normalize to the reference).}
  \item{f.cdf}{a \code{Cdf} object.}
  \item{ref}{a vector of reference intensities.}
  \item{prd.td}{cutoff parameter (details in the bibliographic reference)}
  \item{progress}{A flag to have a dumps throughout the normalization}
}
\value{
  Respectively a \code{\link{Cel.container}} of normalized objects, or
  a vector of normalized intensities, with an attribute "invariant.set"
  holding the indexes of the 'invariant' intensities.
}
\details{
  The set of invariant intensities between \code{data} and
  \code{ref} is found through an iterative
  process (based on the respective ranks the intensities).
  This set of intensities is used to generate a normalization curve by
  smoothing.
}
\author{
  L. Gautier <laurent@cbs.dtu.dk>
  (Thanks to Cheng Li for the discussions about the algorithm.)
}
\references{
  Cheng Li and Wing Hung Wong, Model-based analysis of
  oligonucleotides arrays: model validation, design issues and
  standard error application. Genome Biology 2001, 2(8):research0032.1-0032.11
}
\seealso{
  \code{\link{normalize}} to normalize \code{\link{Cel.container}} objects.
}
%\examples{
%  
%}
\keyword{manip}
