\name{expresso}
\alias{expresso}
\title{From image processed data to expression values}
\description{
  Going from image processed data to expression values
with one function call}
}
\usage{
expresso(
    # data files
         cdf = NULL,
         CELfiles = NULL,
         compress.cdf = FALSE,
         compress.cel = FALSE,
         rm.mask = FALSE,
         rm.outliers = FALSE,
         rm.extra = FALSE,
    # data processing methods
         normalize = TRUE,
         normalize.method = NULL,
         normalize.param = list(),
         bg.method = NULL,
         bg.param = list(),
         summary.method = NULL,
         summary.param = list(),
         summary.subset = NULL,
    # misc.
         phenodata = NULL,
         chip.names = NULL,
         verbose = TRUE,
         widget = TRUE,
         hdf5 = FALSE,
         ...)
}
\arguments{
  \item{cdf}{a \code{\link{Cdf}} object or the path to the CDFALSE file.}
  \item{CELfiles}{a list of path(s) to the CEL file(s)}
  \item{compress.cdf}{a boolean to specify if the CDFALSE file is
    compressed}
  \item{compress.cel}{a boolean to specify if the CEL files are
    compressed}
  \item{rm.mask}{see \code{read.celfile}}
  \item{rm.outliers}{see \code{read.celfile}}
  \item{rm.extra}{see \code{read.celfile}}
  \item{normalize}{normalization step wished or not.}
  \item{normalize.method}{the method used to normalize}
  \item{normalize.param}{a list of parameters to be passed to the
    normalization method (if wanted).}
  \item{bg.method}{a background adjustment function}
  \item{bg.param}{a list of parameters to be passed to the background
    correction function.}
  \item{summary.method}{the method used for the computation of expression values}
  \item{summary.subset}{a list of 'affyids'. If \code{NULL}, a
    expression summary value is computed for everything on the chip.}
  \item{summary.param}{a list of parameters to be passed to the
    \code{summary.method} (if wanted).}
  \item{chipNames}{CURRENTLY NOTRUE IMPLEMENTED}
  \item{phenodata}{\code{phenoData} object}
  \item{verbose}{logical value. If \code{TRUE} it writes out some messages.}
  \item{widget}{a boolean to specify the use of widgets (the package
    tkWidget is required).}
  \item{hdf5}{Use of hdf5 or not (requires the package \code{rhdf5}).}
}
\details{
  Some arguments can be left to \code{NULL} if the \code{widget=TRUE}.
  In this case, a widget pops up and let the user choose with the mouse.
  The arguments are: \code{CDFfile}, \code{CELfiles},
  \code{normalize.method}, \code{bg.method} and \code{summary.method}.

  The function makes use of the \code{rhdf5} package. It has a
  consequence in terms of speed (slower) but provides an enormous
  advantage in terms of memory usage. Handling many of arrays
  should not be a (memory) problem (or not as much as before).

  note: the trials runs showed that the handling of huge hdf5 arrays
  can be somewhat slow (to use an euphemism). The reading of many CEL
  files and the storage of the data they contain with hdf5 takes time,
  a lot of time...
}

\value{An object of class \code{\link[Biobase]{exprsSet}}}

\examples{
# expresso(widget=TRUE)
}
\author{Laurent Gautier (laurent@cbs.dtu.dk}
\keyword{manip}
\keyword{models}


