\name{xy2indices}
\alias{xy2indices}
\alias{indices2xy}
\title{Functions to convert indices to x/y (and reverse)}
\description{
  Functions to convert indices to x/y (and reverse)
}
\usage{
xy2indices(x, y, nc = NULL, cel = NULL, abatch = NULL, cdf = NULL, xy.offset = NULL)
indices2xy(i, nc = NULL, cel = NULL, abatch = NULL, cdf = NULL, xy.offset = NULL)
}
\arguments{
  \item{x}{A numeric vector of \code{X} (column) position(s) for the probes.}
  \item{y}{A numeric vector of \code{Y} (row) position(s) for the probes.}
  \item{i}{A numeric vector of indices in the \code{AffyBatch} for the probes.}
  \item{nc}{total number of columns on the chip. It is usually better to
  specify either the cdf or abatch arguments rather than the number of columns.}
  \item{cel}{a corresponding object of class \code{Cel}. This has been
    deprecated. Use abatch or cdf instead.}
  \item{abatch}{a corresponding object of class
    \code{\link[affy:AffyBatch-class]{AffyBatch}}.}
  \item{cdf}{character - the name of the corresponding cdf package.}
  \item{xy.offset}{an eventual offset for the XY coordinates. See Details.}
}
\details{
  The Affymetrix scanner reads data from a GeneChip by row, and exports
  those data to a CEL file. When we read in the CEL file data to an
  \code{AffyBatch} object, we store data for each GeneChip as a single
  column in a matrix of probe-wise intensity values.

  The CDF files that Affymetrix make available for various GeneChips map
  individual probes to probesets based on their (x,y) coordinates on the
  GeneChip. Note that these coordinates are zero-based, and (x,y) is the
  same as (column, row). In other words, the x coordinate indicates the
  horizontal location of the probe, and the y coordinate indicates the
  vertical location of the probe. By convention, (0,0) is the coordinate
  location for the top left position, and (ncol-1, nrow-1) is the
  coordinate location of the lower right position.

  For most users, the mapping of probes to probeset is handled
  internally by various functions (\code{rma}, \code{espresso}, etc),
  and in general usage it is never necessary for a user to convert probe
  index position in an \code{AffyBatch} to the corresponding (x,y)
  coordinates on the GeneChip. These functions are only useful for those
  who wish to know more about the internal workings of the Affymetrix
  GeneChip.
  
  The parameter \code{xy.offset} is there for compatibility.
  For historical reasons, the xy-coordinates for the features
  on Affymetrix GeneChips were decided to start at 1 (one) rather than 0
  (zero). One can set the offset to 1 or to 0. Unless the you \_really\_
  know what you are doing, it is advisable to let it at the default
  value \code{NULL}. This way the package-wide option \code{xy.offset} is
  always used.
}
\value{
  A vector of indices or a two-columns matrix of Xs and Ys. 
}
\author{L.}
\section{Warning}{Even if one really knows what is going on, playing with
  the parameter \code{xy.offset} could be risky. Changing the package-wide
  option \code{xy.offset} appears much more sane.}
\seealso{\code{\link{indexProbes}}}
\examples{
if (require(affydata)) {
  data(Dilution)
  pm.i <- indexProbes(Dilution, which="pm", genenames="AFFX-BioC-5_at")[[1]]
  mm.i <- indexProbes(Dilution, which="mm", genenames="AFFX-BioC-5_at")[[1]]

  pm.i.xy <- indices2xy(pm.i, abatch = Dilution)
  mm.i.xy <- indices2xy(mm.i, abatch = Dilution)

  ## and back to indices
  i.pm <- xy2indices(pm.i.xy[,1], pm.i.xy[,2], cdf = "hgu95av2cdf")
  i.mm <- xy2indices(mm.i.xy[,1], mm.i.xy[,2], cdf = "hgu95av2cdf")

  identical(pm.i, as.integer(i.pm))
  identical(mm.i, as.integer(i.mm))

  image(Dilution[1], transfo=log2)
  ## plot the pm in red
  plotLocation(pm.i.xy, col="red")
  plotLocation(mm.i.xy, col="blue")
}
}
\keyword{manip}
