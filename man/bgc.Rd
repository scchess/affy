\name{bg.correct}
\alias{bg.correct}
\alias{bg.correct-methods}
\alias{bg.correct.methods}
\alias{bg.correct.none} %took out .pmonly casue rma is pm-only
%\alias{bg.correct.subtractmm}
\alias{bg.correct.rma}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Background Correction}
\description{
Background corrects the probe values in an object of class
\code{\link{AffyBatch}}.
}
\usage{
bg.correct(object, method)

%bg.correct.subtractmm(object)
bg.correct.Affy.chipwide

bg.correct.rma(object)

bg.correct.none(object)}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{object}{An object of class \code{\link{AffyBatch}}.}
  \item{method}{An element of mode \code{character} defining one of the
    methods defined in \code{bg.correct.methods}. }
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

  Details on the internal fucntions used by \code{bg.correct.rma} are in
  \code{\link{bg.adjust}}.
  
}  
\value{
  A \code{\link{AffyBatch}} for which the background has been
  corrected. For some methods(*), only PMs are corrected and the MMs remain
  the same.
  (*) rma.
}
\references{Irizarry et al. (2002)}

\examples{
    data(affybatch.example)
    ##the bg corrected AffyBatch 
    bgc <- bg.correct(affybatch.example,method="rma") 
    ##the same as bg.correct.rma(affybatch.example)
    plot(pm(affybatch.example)[,1],pm(bgc)[,1],log="xy",
	main="PMs before and after background correction")
    #no subtract MM
## chipwide segfaults (somewhere in the C code)
    #bgc <- bg.correct(affybatch.example,method="Affy.chipwide")
    #Index <- pm(bgc)[,1]>0
    #plot(pm(affybatch.example)[Index,1],pm(bgc)[Index,1],log="xy",
	#main="PMs before and after background correction")
 
}
\keyword{manip}

