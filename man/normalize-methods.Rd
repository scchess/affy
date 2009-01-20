\name{normalize-methods}
\title{Normalize Affymetrix Probe Level Data - methods}
\docType{methods}
\alias{normalize.AffyBatch}
\alias{normalize.methods}
\alias{normalize.AffyBatch.methods}
\alias{upDate.normalize.AffyBatch.methods}
\alias{normalize,AffyBatch-method}
\alias{normalize.methods,AffyBatch-method}
\alias{normalize.AffyBatch.methods}
\alias{bgcorrect.methods}
\alias{upDate.bgcorrect.methods}
\alias{pmcorrect.methods}
\alias{upDate.pmcorrect.methods}

\description{
  Method for normalizing Affymetrix Probe Level Data
}
\usage{
normalize.methods(object)
bgcorrect.methods()
upDate.bgcorrect.methods(x)
pmcorrect.methods()
upDate.pmcorrect.methods(x)
}
\arguments{
  \item{object}{An \code{\link[affy:AffyBatch-class]{AffyBatch}}.}
  \item{x}{A character vector that will replace the existing one.}
}
\details{

  If \code{object} is an
  \code{\link[affy:AffyBatch-class]{AffyBatch}} then
  \code{normalize(object)} returns an
  \code{\link[affy:AffyBatch-class]{AffyBatch}} with the intensities
  normalized using the methodology specified by
  \code{getOption("BioC")$affy$normalize.method}. The affy package
  default is \code{quantiles}.

  Other methodologies can be used by specifying them with the
  \code{method} argument. For example to use the invariant set
  methodology described by Li and Wong (2001) one would type:
  \code{normalize(object, method="invariantset")}.

  Further arguments passed by \code{...}, apart from \code{method},  are
  passed along to the 
  function responsible for the methodology defined by the \code{method}
  argument. 
  
  A character vector of \emph{nicknames} for the 
  methodologies available is returned by
  \code{normalize.methods(object))}, where \code{object} 
  is an  \code{\link[affy:AffyBatch-class]{AffyBatch}}, or simply by
  typing \code{normalize.AffyBatch.methods}. If the nickname of a method
  is called "loess", the help page for that 
  specific methodology can
  be accessed by typing \code{?normalize.loess}. 

  For more on the normalization methodologies currently implemented
  please refer to the vignette `Custom Processing Methods'.

  To add your own normalization procedures please
  refer to the customMethods vignette.


  The functions: \code{bgcorrect.methods},
\code{pmcorrect.methods}, provide access to internal vectors
listing the corresponding capabilities.

}
\seealso{
  \code{\link{AffyBatch-class}}, \code{\link{normalize}}.
}
\examples{
if (require(affydata)) {
  data(Dilution)
  normalize.methods(Dilution)
  generateExprSet.methods()
  bgcorrect.methods()
  pmcorrect.methods()
}
}
\keyword{manip}
