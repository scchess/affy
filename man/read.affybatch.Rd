\name{read.affybatch}
\alias{read.affybatch}

\title{ Read CEL files into an AffyBatch}
\description{
  Read CEL files into an Affybatch
}
\usage{
read.affybatch(..., filenames = character(0),
               phenoData = new("phenoData"), description = new("miame"),
               notes = "",
               compress = getOption("BioC")$affy$compress.cel,
               rm.mask = FALSE, rm.outliers = FALSE, rm.extra = FALSE,
               hdf5 = FALSE, hdf5FilePath = NULL,
               verbose = FALSE)
}
\arguments{
  \item{\dots}{ file names }
  \item{filenames}{ file names }
  \item{phenoData}{ a \code{phenoData} object}
  \item{description}{ a \code{miame} object }
  \item{notes}{ notes }
  \item{compress}{ are the CEL files compressed ?}
  \item{rm.mask}{ should the spots marked as 'MASKS' set to \code{NA} ? }
  \item{rm.outliers}{ should the spots marked as 'OUTLIERS' set to \code{NA} }
  \item{rm.extra}{ if \code{TRUE}, overrides what is in \code{rm.mask} and \code{rm.oultiers} }
  \item{hdf5}{ use of hdf5 ? }
  \item{hdf5FilePath}{ a filename to use with hdf5 }
  \item{verbose}{ verbosity flag }
}
\details{
  ...stay tuned...
}
\value{
  An \code{AffyBatch} object.
}

\seealso{\code{\link{AffyBatch}}}

\examples{

}
\keyword{manip}
