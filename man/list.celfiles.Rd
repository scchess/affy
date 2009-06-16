\name{list.celfiles}
\alias{list.celfiles}
\title{List the Cel Files in a Directory/Folder}
\description{
This function produces a vector containing the names of files in the
named directory/folder ending in .cel or .CEL. 
}
\usage{
list.celfiles(...)
}
\arguments{
  \item{\dots}{arguments to  pass along to
    \code{\link[base]{list.files}}}
}
\value{
A character vector of file names.
}
\seealso{list.files}
\examples{
list.celfiles()
}
\keyword{character}
