\name{cleancdfname}
\alias{cleancdfname}
\alias{mapCdfName}
\title{Clean Affymetrix's CDF name}
\description{
This function converts Affymetrix's names for CDF files to the names
used in the annotation package and in all Bioconductor.
}
\usage{
cleancdfname(cdfname, addcdf = TRUE)
}
\arguments{
  \item{cdfname}{A \code{character} denoting Affymetrix'x CDF file name }
  \item{addcdf}{A \code{logical}. If \code{TRUE} it adds the string "cdf"
    at the end of the cleaned CDF name. This is used to name the
    \code{cdfenvs} packages.}
}
\details{
  This function takes a CDF filename obtained from an Affymetrix file
  (from a CEL file for example) and convert it to a convention of ours:
  all small caps and only alphanumeric characters. The details of the rule
  can be seen in the code.
  We observed exceptions that made us create a set of special cases for
  mapping CEL to CDF. The object \code{mapCdfName} holds information
  about these cases. It is a \code{data.frame} of three elements: the 
  first is the name as found in the CDF file, the second the name in the
  CEL file and the third the name in Bioconductor. \code{mapCdfName} can
  be loaded using \code{data(mapCdfName)}.
}
\value{
A \code{character}
}
\examples{
cdf.tags <- c("HG_U95Av2", "HG-133A")
for (i in cdf.tags)
  cat(i, "becomes", cleancdfname(i), "\n")
}
\keyword{character}
