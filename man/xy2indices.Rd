\name{xy2indices}
\alias{xy2indices}
\alias{indices2xy}
\title{ Functions to convert indices to x/y (and reverse)}
\description{
  Functions to convert indices to x/y (and reverse)
}
\usage{
xy2indices(x, y, nr = NULL, cel = NULL, abatch = NULL)
indices2xy(i, nr = NULL, cel = NULL, abatch = NULL)
}
\arguments{
  \item{x}{ \code{X} position for the probes }
  \item{y}{ \code{Y} position for the probes }
  \item{i}{ indices in the \code{AffyBatch} for the probes }
  \item{nr}{ total number of \code{Xs} on the chip }
  \item{cel}{ a corresponding object of class \code{\link[affy:Cel-class]{Cel}} }
  \item{abatch}{ a corresponding object of class \code{\link[affy:AffyBatch-class]{AffyBatch}} }
}
\details{
  The probes intensities for given probe set ids are extracted from an
  \code{AffyBatch} object using the indices stored in the corresponding
  \code{cdfenv}.
}
\value{
  A vector of indices or a two-columns matrix of Xs and Ys. 
}
\author{ L. }
\warnings{  The X/Y values stored in {\it Affymetrix} CDF and probe sequence files
  start their indexing at zero. One should be careful and add one
  to these values\ldots a unified scheme with \code{oligochips} might
  come in the future (release 1.3), and make the (dangerous) juggling
  with zero-offset/one-offset unnecessary\ldots~~ }
\seealso{\code{\link{locate.name}}, \code{\link{indexProbes}}}
\examples{
data(affybatch.example)

pm.i <- indexProbes(affybatch.example, which="pm", genenames="AFFX-BioC-5_at")[[1]]
mm.i <- indexProbes(affybatch.example, which="mm", genenames="AFFX-BioC-5_at")[[1]]

pm.i.xy <- indices2xy(pm.i, abatch = affybatch.example)
mm.i.xy <- indices2xy(mm.i, abatch = affybatch.example)

image(affybatch.example[[1]], transfo=log2)
## plot the pm in red
plotLocation(pm.i.xy, col="red")
plotLocation(mm.i.xy, col="blue")


}
\keyword{ manip }
