\name{rma}
\alias{rma}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{RMA}
\description{
This function computes the RMA measure.}
}
\usage{
rma(object, subset = NULL, phenodata = NULL, annotation = NULL, description = NULL, notes = NULL, verbose = TRUE, ...)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{object}{Plob}
 \item{subset}{a boolean or indexing vector indicating a subset of gene
    names (see \code{geneNames}).}
  \item{phenodata}{a \code{phenoData} object to associate to the
    \code{exprSet} returned.}
  \item{annotation}{annotation to associate to the \code{exprSet} returned.}
  \item{description}{description to associate to the \code{exprSet} returned.}
  \item{notes}{notes to associate to the \code{exprSet} returned.}
  \item{verbose}{logical value. If \code{TRUE} it writes out some messages.}
  \item{...}{parameters passed to the \code{summary.stat} function given.}
}
\details{
This function is like \code{\link{express}} but it is hard wired to
compute RMA}
\value{An object of class \code{\link[Biobase]{exprsSet}}}

\author{Rafael A. Irizarry, Ben Bolstad, Leslie Cope}
\seealso{\code{\link{express}}}
  
\keyword{manip}% at least one, from doc/KEYWORDS
\keyword{models}% __ONLY ONE__ keyword per line
