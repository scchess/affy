\name{AffyBatch-class}
\docType{class}
\alias{AffyBatch-class}
\alias{AffyBatch}
\title{Class "AffyBatch"}
\description{ A class to bundle data from several experiments }
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("AffyBatch", ...)}.
    ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{intensity}:}{Object of class \code{"array"}. Contains
      the probe intensities}
    \item{\code{cdfName}:}{Object of class \code{"character"}. The name
      of the related CDF file.}
    \item{\code{chipNames}:}{Object of class \code{"character"}. The
      names of the chips in the batch. }
    \item{\code{phenoData}:}{Object of class \code{"phenoData"} ~~ }
    \item{\code{nrow}:}{Object of class \code{"numeric"} The number of
      rows for each chip in the batch }
    \item{\code{ncol}:}{Object of class \code{"numeric"} The number of
      columns for each chip in the batch }
    \item{\code{nexp}:}{Object of class \code{"numeric"}. The number of
      experiments in the batch.}
    \item{\code{annotation}:}{Object of class \code{"character"} ~~ }
    \item{\code{description}:}{Object of class \code{"characterORmiame"} ~~ }
    \item{\code{notes}:}{Object of class \code{"character"} ~~ }
    \item{\code{history}:}{Object of class \code{"list"} ~~ }
  }
}

\section{Methods}{
  \describe{
    \item{[[<-}{\code{signature(x = "AffyBatch")}: ... }
    \item{[[}{\code{signature(x = "AffyBatch")}: ... }
    \item{chipNames}{\code{signature(object = "AffyBatch")}: ... }
    \item{history}{\code{signature(max.show = "AffyBatch", reverse = "missing")}: ... }
    \item{history<-}{\code{signature(object = "AffyBatch")}: ... }
    \item{intensity<-}{\code{signature(object = "AffyBatch")}: ... }
    \item{intensity}{\code{signature(object = "AffyBatch")}: ... }
    \item{length}{\code{signature(x = "AffyBatch")}: ... }
    \item{normalize}{\code{signature(object = "AffyBatch")}: ... }
    \item{normalize.methods}{\code{signature(object = "AffyBatch")}: ... }
    \item{show}{\code{signature(object = "AffyBatch")}: ... }
    \item{chipNames<-}{\code{signature(object = "AffyBatch")}: ... }
    \item{computeExprSet}{\code{signature(x = "AffyBatch", bg.method = "character", summary.method = "character")}: ... }
    \item{geneNames}{\code{signature(object = "AffyBatch")}: ... }
    \item{getCdfInfo}{\code{signature(object = "AffyBatch")}: ... }
    \item{indexProbes}{\code{signature(object = "AffyBatch", which = "character")}: ... }
    \item{mm}{\code{signature(object = "AffyBatch")}: ... }
    \item{mmindex}{\code{signature(object = "AffyBatch")}: ... }
    \item{pm}{\code{signature(object = "AffyBatch")}: ... }
    \item{pmindex}{\code{signature(object = "AffyBatch")}: ... }
    \item{ppset}{\code{signature(object = "AffyBatch")}: ... }
    \item{probeNames}{\code{signature(object = "AffyBatch")}: ... }
    \item{probes}{\code{signature(object = "AffyBatch")}: ... }
  }
}
\seealso{ supplementary methods for \code{\link{AffyBatch-class}}, are
  \code{\link{pairs.AffyBatch}}, \code{\link{merge.AffyBatch}} and
  \code{\link{expresso.AffyBatch}} }

\examples{
data(affybatch.example)
affybatch.example
}
\keyword{classes}
