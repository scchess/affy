\name{bg.correct}
\alias{bg.correct}
\alias{bg.correct-methods}
\alias{bg.correct.methods}
\alias{bgcorrect.methods}
\alias{bg.correct.none} %took out .pmonly casue rma is pm-only
%\alias{bg.correct.subtractmm}
\alias{bg.correct.rma}
\alias{bg.correct.rma2}
\alias{bg.correct.mas}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Background Correction}
\description{
Background corrects the probe values in an object of class
\code{\link[affy:AffyBatch-class]{AffyBatch}}.
}
\usage{
bg.correct(object, method, ...)

%bg.correct.subtractmm(object)
bg.correct.rma(object,...)
bg.correct.mas(object, griddim)

bg.correct.none(object, ...)}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{object}{An object of class \code{\link[affy:AffyBatch-class]{AffyBatch}}.}
  \item{method}{An element of mode \code{character} defining one of the
    methods defined in \code{bg.correct.methods}. }
  \item{griddim}{grid dimension used for background estimate. The array
    is divided into griddm equal parts. Default is 16.}
  \item{...}{arguments to pass along to the engine function.}
}
\details{
  The name of the method to apply must be double-quoted.
  Methods provided with the package are currently:
  \begin{itemize}
  \item bg.correct.none: return the object.
  %  \item bg.correct.subtractmm: subtract the MMs from the PMs.
  \item bg.correct.chipwide: noise correction as described in the `sadd
  white paper' from Affymetrix.
  \item bg.correct.rma: the model based correction used by the rma
  expression measure.
  \end{itemize}
  They are listed in the variable bg.correct.methods. The user must
  supply the word after "bg.correct", i.e none, subtractmm, rma, etc...

  More details are available in the vignette.

  \code{bg.correct.rma2} gives access to the older background method
  that was used by the \code{\link{rma}} function in the version 1.1
  bioconductor releases. For the 1.2 release \code{\link{rma}} performs
  the same background adjustment as \code{bg.correct.rma}. Unless you
  wish to duplicate older results you should use \code{bg.correct.rma}
  
  Details on the internal fucntions used by \code{bg.correct.rma} are in
  \code{\link{bg.adjust}}.
  
}  
\value{
  A \code{\link[affy:AffyBatch-class]{AffyBatch}} for which the background has been
  corrected. For some methods (rma), only PMs are corrected and the MMs remain
  the same.
}
\references{Irizarry et al. (2002)}

\examples{
    data(affybatch.example)
    ##the bg corrected AffyBatch 
    bgc <- bg.correct(affybatch.example,method="rma") 
    ##the same as bg.correct.rma(affybatch.example)
    plot(pm(affybatch.example)[,1],pm(bgc)[,1],log="xy",
	main="PMs before and after background correction")

    #the mas 5.0 background
    #bgc <- bg.correct(affybatch.example,method="mas")
    #Index <- pm(bgc)[,1]>0
    #plot(pm(affybatch.example)[Index,1],pm(bgc)[Index,1],log="xy",
	#main="PMs before and after background correction")
 
}
\keyword{manip}

