\name{whatcdf}
\alias{whatcdf}

\title{Find which CDFALSE corresponds}
\description{
  Find which kind of CDFALSE corresponds to a CEL file.
}
\usage{
whatcdf(filename, compress = FALSE)
}
\arguments{
  \item{filename}{a '.CEL' file name}
  \item{compress}{boolean (file compressed or not}
}
\details{
  Information concerning the corresponding CDFALSE file seems to be found in
  CEL files. This allows us to try to link CDFALSE information automatically.
}
\value{
  a \code{character} with the name of the CDF
}
\seealso{\code{getInfoInAffyFile}, \code{read.celfile}}

\keyword{manip}

