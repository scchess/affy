\name{normalize.invariantset}
\alias{normalize.invariantset}
\title{Scale microarray data}
\description{
  normalizes two microarray data to each others.
}
\usage{
normalize.invariantset(data, ref, prd.td=c(0.003,0.007))
}
\arguments{
  \item{data}{a vector of intensities on a chip (to normalize to the reference).}
  \item{ref}{a vector of reference intensities.}
  \item{prd.td}{A range for threshold values.}
}
\value{
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
  Laurent\cr
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
%  ## TO DO...
%}
\keyword{}


