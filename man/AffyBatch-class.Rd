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

%   abatch.n <- normalize(affybatch.example, CDF.example, method="quantiles")
%   par(mfrow=c(2,1))
%   plot(intensity(affybatch.example[[1]]), intensity(affybatch.example[[2]]), col="black",
%        xlab=sampleNames(affybatch.example)[1], ylab=sampleNames(affybatch.example[[2]]),log="xy")
%   plot(intensity(abatch.n[[1]]), intensity(abatch.n[[2]]), col="red",
%        xlab=sampleNames(affybatch.example[[1]]), ylab=sampleNames(affybatch.example[[2]]),log="xy")

%   norm.method <- c("constant", "invariantset" , "quantiles")

%   par(mfrow=c(3,3))

%   ## loop over the methods
%   for (i in 1:3) {

%     cat("--->",norm.method[i],"\n")  
%     abatch.n <- normalize(affybatch.example, method=norm.method[i])

%     ## plot the raw and normalized intensities
%     plot(c(intensity(affybatch.example[[1]]), intensity(abatch.n[[1]])),
%          c(intensity(affybatch.example[[2]]), intensity(abatch.n[[2]])),
%          xlab=sampleNames(affybatch.example[[1]]), ylab=sampleNames(affybatch.example[[2]]),
%          type="n", main=paste("method:", norm.method[i]))
%     points(intensity(affybatch.example[[1]]), intensity(affybatch.example[[2]]), col="black")
%     points(intensity(abatch.n[[1]]), intensity(abatch.n[[2]]), col="red")
%     ## indentity line
%     abline(0,1)

%     ## kernel density estimates of the intensities
%     md1 <- density(intensity(affybatch.example[[1]]), na.rm=TRUE)
%     md2 <- density(intensity(affybatch.example[[2]]), na.rm=TRUE)
%     n.md1 <- density(intensity(abatch.n[[1]]), na.rm=TRUE)
%     n.md2 <- density(intensity(abatch.n[[2]]), na.rm=TRUE)

%     ## plot the raw an normalized intensities (log scale)
%     plot(range(md1$x, md2$x, n.md1$x, n.md2$x),
%          range(md1$y, md2$y, n.md1$y, n.md2$y),
%          xlab="intensity (log scale)", ylab="density estimate",
%          type="n", log="x")
%     points(md1, type="l", lty="solid", col="black")
%     points(n.md1, type="l", lty="solid", col="red")
%     points(md2, type="l", lty=5, col="black")
%     points(n.md2, type="l", lty=5, col="red")

%     plot(c(intensity(affybatch.example[[1]]), intensity(abatch.n[[1]])),
%          c(intensity(affybatch.example[[2]]), intensity(abatch.n[[2]])),
%          xlab=sampleNames(affybatch.example[[1]]), ylab=sampleNames(affybatch.example[[2]]), log="xy",
%          type="n", main=paste(norm.method[i]," - axis in log scale"))
%     points(intensity(affybatch.example[[1]]), intensity(affybatch.example[[2]]), col="black")
%     points(intensity(abatch.n[[1]]), intensity(abatch.n[[2]]), col="red")
%   }

   panel.idline <- function(x, y, ...) {
     abline(0,1, col="gray")
     points(x, y, ...)
   }
  pairs(affybatch.example, panel=panel.idline, pch=".")

}
\keyword{classes}
