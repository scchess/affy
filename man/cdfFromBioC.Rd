\name{cdfFromBioC}
\alias{cdfFromBioC}
\alias{cdfFromLibPath}
\alias{cdfFromEnvironment}
\title{Functions to obtain CDF files}
\description{
  A set of functions to obtain CDF files from various locations.
}
\usage{
cdfFromBioC(cdfname, lib = .libPaths()[1], verbose = TRUE)
cdfFromLibPath(cdfname, lib = NULL, verbose=TRUE)
cdfFromEnvironment(cdfname, where, verbose=TRUE)
}
\arguments{
  \item{cdfname}{name of the CDF.}
  \item{lib}{install directory for the CDF package.}
  \item{where}{environment to search.}
  \item{verbose}{logical controlling extra output.}  
}
\details{
  These functions all take a requested CDF environment name and will
  attempt to locate that environment in the appropriate location (a
  package's data directory, as a CDF package in the .libPaths(), from a
  loaded environment or on the Bioconductor website.  If the environment
  can not be found, it will return a list of the methods tried that failed.
}
\value{
  The CDF environment or a list detailing the failed locations.
}
\author{Jeff Gentry}
\keyword{utilities}
