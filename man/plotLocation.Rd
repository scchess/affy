\name{plotLocation}
\alias{plotLocation}
\title{Plot a location on a cel image}
\description{
  Plots a location on a previously plotted cel image. This can be used
  to locate the physical location of probes on the array.
}
\usage{
plotLocation(x, cdf=NULL, nrow=NULL, ncol=NULL,
              col="green", pch=22, ...)
}
\arguments{
  \item{x}{a location. It can be obtained by the function
    get.location, or made elsewhere (basically a location is nrows and
    two columns array. The first column corresponds to the x positions
    and the second columns corresponds to the y positions of n elements
    to locate)}
  \item{cdf}{A \code{\link{Cdf}} object. It is used only to scale the
    locations to 
  the existing image.}
  \item{nrow}{The number of rows. Ignored if a \code{\link{Cdf}} object
    is specied.} 
  \item{ncol}{The number of columns. Ignored if a \code{\link{Cdf}}
    object is specied.} 
  \item{col}{colors for the plot}
  \item{pch}{plotting type (see function \code{plot})}
  \item{...}{Other parameters passed to the function \code{points}}
  }
}
\author{
  Laurent Gautier
}
\examples{
  ## loading data
  data(listcel)
  data(CDF.example)

  ## image of the celfile
  image(listcel[[1]])
  
  ## getting the location for the one called "A28102_at"
  l <- locate.name("A28102_at",CDF.example)

  ## plot
  plotLocation(l, CDF.example)

  ## the same without passing the cdffile object as a parameter
  ## (it may prove usefull to save memory)
  f.row <- nrow(CDF.example@name)
  f.col <- ncol(CDF.example@name)
  plotLocation(l,nrow=f.row,ncol=f.col)
}
\keyword{aplot}



