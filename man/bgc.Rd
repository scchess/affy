\name{bgc}
\alias{bgc}
\alias{bg.correct}
\alias{bg.correct-methods}
\alias{bg.correct.methods}
\alias{bg.correct.pmonly}
\alias{bg.correct.subtractmm}
%\alias{bg.correct.adjust}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Background Correction}
\description{
Background corrects the PM values in an object of class
\code{\link{Plob}} or \code{\link{PPSet.container}} using an 
arbitraty function.
}
\usage{
bgc(object, bg = bg.adjust)

bg.correct(object, method)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{object}{An object of class \code{\link{Plob}}.}
  \item{bg}{a function that determines how background correction is
    done. The default is \code{\link{bg.adjust}}.}
  \item{method}{An element of mode \code{character} designating the name
    of the function.
    %DEBUG: still a name ?
  }
}
\details{
  The function that performs the actual background correction is an
  arbitrary function. The implementation differs slightly between
  \code{Plob} and \code{PPSet.container}, but the principle is similar.
  In both case the two first arguments must be PMs values and MMs
  values. In the case of \code{Plob}, they are vectors of all the PMs
  and MMs on one arry, while for \code{PPSet.container} they are
  matrices (one row per probe and one column per array.

  details specific to \code{Plob}:
  The default function is \code{\link{bg.correct.adjust}}.
  The methods provided with the package are:
  \begin{itemize}
  \item \code{subtractmm}: subtracts the MMs
  \item \code{background.adjust}: please refer to the corresponding help file.
  \end{itemize}

  details specific to \code{PPSet.container}:
  The name of the function to apply must be double-quoted.
  Methods provided with the package are currently:
  \begin{itemize}
  \item bg.correct.pmonly: only take the PMs and forget about the MMs
  \item bg.correct.subtractmm: subtract the MMs to the PMs
  \item bg.correct.adjust: worth lengthy explanations. Please
  refer to the corresponding help file.
  \end{itemize}
  They are listed in the variable bg.correct.methods

  (sorry for the mess, we are under development and trying to get/preserve a life too
  ;) )
  
}  
\value{
  A \code{\link{Plob}}
  (or a matrix of `PM corrected values' for PPSet.container)
}
\references{Irizarry et al. (2002)}
\author{
  Rafael A. Irizarry\cr
  and El Gringo (Laurent)}

\seealso{ \code{\link{bg.adjust}}}
\examples{
     data(Dilution)
     tmp <- bgc(Dilution)
     plot(pm(Dilution)[,1],pm(tmp)[,1],log="xy",
	  main="PMs before and after background correction")
}
\keyword{manip}

