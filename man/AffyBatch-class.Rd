\name{AffyBatch-class}
\docType{class}
\alias{AffyBatch-class}
\alias{AffyBatch}
\alias{probes}
\alias{computeExprSet}
\alias{geneNames}
\alias{getCdfInfo}
\alias{image}
\alias{indexProbes}
\alias{intensity<-}
\alias{intensity}
\alias{pmindex}
\alias{mmindex}
\alias{probeset}

\title{Class AffyBatch}
\description{This is a class representation for Affymetrix GeneChip probe
  level data. The main component are the intensities from multiple
  arrays of the same \code{CDF} type. In extends \code{\link[Biobase]{exprSet-class}}.}
\section{Objects from the Class}{
  Objects can be created using the function \code{\link{read.affybatch}}
  or the wrapper \code{\link{AffyBatch}}.}
\section{Slots}{
  \describe{
    \item{\code{cdfName}:}{Object of class \code{"character"} representing the name of
      \code{CDF} file associated with the arrays in the \code{AffyBatch}.}
    \item{\code{nrow}:}{Object of class \code{"numeric"} representing
      the physical number of rows in the arrays.}
    \item{\code{ncol}:}{Object of class \code{"numeric"} representing
      the physical number of columns in the arrays.}
    \item{\code{exprs}:}{Object of class "matrix". The observed
      intensities. This is 
     a matrix with columns representing arrays and rows
     representing cel intensities. }
    \item{\code{se.exprs}:}{Object of class "matrix". This is
     a matrix of the same dimensions as \code{exprs} which may contains
     standard error estimates for the estimated intensity levels. }
    \item{\code{phenoData}:}{Object of class "phenoData" This is an
      instance of class \code{phenoData} 
    containing the patient (or case) level data. The columns of the
    pData slot of this entity represent
    variables and the rows represent patients or cases. }
    \item{\code{annotation}}{A character string identifying the
      annotation that may be used for the \code{exprSet} instance.}
    \item{\code{description}:}{Object of class "MIAME". For
      compatibility with previous version of this class description can
      also be a "character". The class \code{characterOrMIAME} has been
      defined just for this.} 
    \item{\code{notes}:}{Object of class "character" Vector of
      explanatory text }
  }
}
\section{Extends}{
Class \code{"exprSet"}, directly.
}
\section{Methods}{
  \describe{
    \item{pm<-}{\code{signature(object = "AffyBatch")}: replaces the
      perfect match intensities }
    \item{pm}{\code{signature(object = "AffyBatch")}: extracts the pm intensities. }
    \item{mm<-}{\code{signature(object = "AffyBatch")}: replaces the
      mismatch intensities.}
    \item{mm}{\code{signature(object = "AffyBatch")}: extracts the mm intensities. }
    \item{probes}{\code{signature(object = "AffyBatch", which)}: extract the
      perfect match or mismatch probe intensities. Uses which can be "pm"
    and "mm".}
    \item{[<-}{\code{signature(x = "AffyBatch")}: replaces subsets. }
    \item{[}{\code{signature(x = "AffyBatch")}: subsets by array.}
    \item{[[}{\code{signature(x = "AffyBatch")}: extracts each array as
      an \code{\link{Cel}} instance.}
    \item{boxplot}{\code{signature(x = "AffyBatch")}:  creates a
      \code{\link{boxplot}s of log base 2 intensities. }
    \item{computeExprSet}{\code{signature(x = "AffyBatch",
	summary.method = "character")}: For each probe set computes an
      expression value using \code{summary.method}. }
    \item{geneNames}{\code{signature(object = "AffyBatch")}: return the
      probe set names also referred to as the Affymetrix IDs. }
    \item{getCdfInfo}{\code{signature(object = "AffyBatch")}: retrieve
      the environment that defines the location of probes by probe set. }
    \item{image}{\code{signature(x = "AffyBatch")}: creates an image for
      each sample.}
    \item{indexProbes}{\code{signature(object = "AffyBatch", which =
	"character")}: returns a list with locations of the probes in
      each probe  set. The list names defines the probe set
      names. \code{which} can be "pm", "mm", or "both". If "both" then
      perfect match locations are given followed by mismatch locations.} 
    \item{intensity<-}{\code{signature(object = "AffyBatch")}: a
      replacement method for the \code{exprs} slot, i.e. the intensities. }
    \item{intensity}{\code{signature(object = "AffyBatch")}: extract the
      \code{exprs} slot, i.e. the intensities. }
    \item{length}{\code{signature(x = "AffyBatch")}: returns the number
      of samples. }
    \item{pmindex}{\code{signature(object = "AffyBatch")}: return the
      location of perfect matches in the intensity matrix.}
    \item{mmindex}{\code{signature(object = "AffyBatch")}: return the
      location of the mismatch intensities. }
    \item{ncol}{\code{signature(x = "AffyBatch")}: An accessor function
      for \code{ncol}.}
    \item{nrow}{\code{signature(x = "AffyBatch")}: an accessor function
      for \code{nrow}}.}
    \item{normalize}{\code{signature(object = "AffyBatch")}: a method to
    \code{\link{normalize}}.}
    \item{normalize.methods}{\code{signature(object = "AffyBatch")}:
      returns the normalization methods defined for this class. See \code{\link{normalize}}.}
    \item{probeNames}{\code{signature(object = "AffyBatch")}: returns
      the probe set associated with each row of the intensity matrix. }
    \item{probeset}{\code{signature(object = "AffyBatch",genenames=NULL,
	locations=NULL)}: Extracts \code{\link{ProbeSet}} objects related to
      the probe sets given in genenames. If an alternative set of
      locations defining pms and mms a list with those locations should
      be passed via the \code{locations} argument.}
    \item{bg.correct}{\code{signature(object = "AffyBatch",
	method="character")}} applies background correction methods
    defined by method.
  }
}
\note{This class is better described in the vignette.}
\seealso{  \code{\link{Cel}}, \code{\link{Cdf}}, 
  \code{\link{exprSet-class}}} 

\examples{
data(affybatch.example)
print(affybatch.example)
pm(affybatch.example)[1:5,]
mm(affybatch.example)[1:5,]
pmindex(affybatch.example,"A28102_at")

data(Dilution)
description(Dilution) ##we can also use the methods of exprSet
sampleNames(Dilution)
abstract(Dilution)
}

\keyword{classes}
