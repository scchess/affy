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
  \item{x}{\code{X} position for the probes.}
  \item{y}{\code{Y} position for the probes.}
  \item{i}{indices in the \code{AffyBatch} for the probes.}
  \item{nc}{total number of columns on the chip.}
  \item{cel}{a corresponding object of class \code{Cel}. This has been
    deprecated. Use abatch or cdf instead.}
  \item{abatch}{a corresponding object of class
    \code{\link[affy:AffyBatch-class]{AffyBatch}}.}
  \item{cdf}{character - the name of the corresponding cdf package.}
  \item{xy.offset}{an eventual offset for the XY coordinates. See Details.}
}
\details{
  The probes intensities for given probe set ids are extracted from an
  \code{AffyBatch} object using the indices stored in the corresponding
  \code{cdfenv}.

  The parameter \code{xy.offset} is there for compatibility.
  For historical reasons, the xy-coordinates for the features
  on Affymetrix chips were decided to start at 1 (one) rather than 0
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

  identical(pm.i, xy2indices(pm.i.xy[,1], pm.i.xy[,2], cdf = "hgu95av2cdf"))
  identical(mm.i, xy2indices(mm.i.xy[,1], mm.i.xy[,2], cdf = "hgu95av2cdf"))

  image(Dilution[1], transfo=log2)
  ## plot the pm in red
  plotLocation(pm.i.xy, col="red")
  plotLocation(mm.i.xy, col="blue")
}
}
\keyword{manip}
