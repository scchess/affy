\name{expresso}
\alias{expresso}

%- Also NEED an `\alias' for EACH other topic documented here.
\title{ From raw probe intensities to expression values }
\description{
  Goes from raw probe intensities to expression values
}
\usage{
expresso(
        afbatch,
    # background correction
        bg.correct = TRUE,
        bgcorrect.method = NULL,
        bgcorrect.param = list(),
    # normalize
        normalize = TRUE,
        normalize.method = NULL,
        normalize.param = list(),
    # pm correction
        pmcorrect.method = NULL,
        pmcorrect.param = list(),
    # expression values
        summary.method = NULL,
        summary.param = list(),
        summary.subset = NULL,
    # misc.
        verbose = TRUE,
        warnings = TRUE,
        widget = FALSE)
}

\arguments{
  \item{afbatch}{ An \code{\link{AffyBatch}} object }
  \item{bg.correct}{a boolean to express whether background correction
    is wanted or not.}
  \item{bgcorrect.method}{the name of the background adjustment method}
  \item{bgcorrect.param}{a list of parameters for bgcorrect.method (if needed/wanted)}
  \item{normalize}{ normalization step wished or not.}
  \item{normalize.method}{the normalization method to use}
  \item{normalize.param}{a list of parameters to be passed to the
    normalization method (if wanted).}
  \item{pmcorrect.method}{the name of the PM adjustement method}
  \item{pmcorrect.param}{a list of parameters for pmcorrect.method (if needed/wanted)}
  \item{summary.method}{the method used for the computation of expression values}
  \item{summary.param}{ a list of parameters to be passed to the
    \code{summary.method} (if wanted).}
  \item{summary.subset}{a list of 'affyids'. If \code{NULL}, a
    expression summary value is computed for everything on the chip.}
  \item{verbose}{logical value. If \code{TRUE} it writes out some
    messages.}
  \item{warnings}{warning when something goes wrong}
  \item{widget}{a boolean to specify the use of widgets (the package
    tkWidget is required).}
}
\details{
    Some arguments can be left to \code{NULL} if the \code{widget=TRUE}.
  In this case, a widget pops up and let the user choose with the mouse.
  The arguments are: \code{AffyBatch},
  \code{bgcorrect.method}, \code{normalize.method}, \code{pmcorrect.method} and \code{summary.method}.

  For the mas 5.0 and 4.0 methods ones need to normalize after obtaining
  expression. The function \code{\link{affy.scalevalue.exprSet}} does this.
}
\value{
  An object of class \code{\link[Biobase]{exprsSet}}
}
%\seealso{\code{\link{}}, \code{\link{}}, \code{\link{}}}
\examples{
data(affybatch.example)

eset <- expresso(affybatch.example, bgcorrect.method="rma",
                 normalize.method="constant",pmcorrect.method="pmonly",
                 summary.method="avgdiff", warnings=FALSE)
}
\keyword{ manip }
