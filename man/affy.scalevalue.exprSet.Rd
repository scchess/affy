%%%THIS FUNCTION SHOULD BE In BIOBASE or SOMEWHERE ELSE... NEXT RELEASE
\name{affy.scalevalue.exprSet}
\alias{affy.scalevalue.exprSet}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Scale normalization for expreSets}
\description{
  Normalizes expression values using the method described in the
  Affymetrix user manual.
}
\usage{
affy.scalevalue.exprSet(eset, sc = 500, analysis="absolute")
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{eset}{An \code{\link{exprSet}} object.}
  \item{sc}{Value at which all arrays will be scaled to.}
  \item{analysis}{should we do absolute or comparison analysis, although
  "comparison" is still not implemented.}
}
\details{
  See the affymetrix user manual for details.
}
\value{
  A normalize \code{\link[Biobase]{exprSet}}}
\keyword{manip}% at least one, from doc/KEYWORDS

