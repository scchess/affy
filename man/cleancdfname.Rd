\name{cleancdfname}
\alias{cleancdfname}
\alias{mapCdfName}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Clean Affymetrix's CDF name}
\description{
This function converts Affymetrix's names for CDF files to the names
used in the annotation pacakge and in all Bioconductor.
}
\usage{
cleancdfname(cdfname, addcdf = TRUE)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{cdfname}{A \code{character} denoting Affymetrix'x CDF file name }
  \item{addcdf}{A \code{logical}. If \code{TRUE} it adds the string
    "cdf" at the end of the cleaned CDF name. This is used to name the
    location packages.}
}
\details{
  This function takes the CDF filenames obtained from the CEL files and
  converts them to the convention used in Bioconductor. All small caps
  and only alphanumeric characters. There is a simple rule, see code,
  but there are also some exception which are defined in
  \code{mapCdfName.tab} in the \code{data} directory.
}
\value{
A \code{character}
}

\examples{
cat("HG_U95Av2 is",cleancdfname("HG_U95Av2"),"\n")
cat("HG-133A is",cleancdfname("HG-133A"),"\n")
}
\keyword{character}% at least one, from doc/KEYWORDS

