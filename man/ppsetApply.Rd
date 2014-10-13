\name{ppsetApply}
\alias{ppsetApply}
\alias{ppset.ttest}
\title{ Apply a function over the ProbeSets in an AffyBatch }
\description{
  Apply a function over the ProbeSets in an AffyBatch
}
\usage{
ppsetApply(abatch, FUN, genenames = NULL, ...)

ppset.ttest(ppset, covariate, pmcorrect.fun = pmcorrect.pmonly, ...)
}
\arguments{
  \item{abatch}{an object inheriting from \code{AffyBatch}.}
  \item{ppset}{an object of class \code{ProbeSet}.}
  \item{covariate}{the name a covariate in the slot \code{phenoData}.}
  \item{pmcorrect.fun}{a function to correct PM intensities.}
  \item{FUN}{a function working on a \code{ProbeSet}.}
  \item{genenames}{a list of Affymetrix probesets ids to work with. All
  probe set ids used when \code{NULL}.}
  \item{\dots}{optional parameters to the function \code{FUN}.}
}
\value{
  Returns a \code{list} of objects, or values, as returned by the
  function \code{FUN} for each \code{ProbeSet} it processes.
}
\author{Laurent Gautier <laurent@cbs.dtu.dk>}
\seealso{\code{\link[affy]{ProbeSet-class}} }
\examples{
ppset.ttest <- function(ppset, covariate, pmcorrect.fun = pmcorrect.pmonly, ...) {
  probes <- do.call("pmcorrect.fun", list(ppset))
  my.ttest <- function(x) {
    y <- split(x, get(covariate))
    t.test(y[[1]], y[[2]])$p.value
  }
  r <- apply(probes, 1, my.ttest)
  return(r)
}
##this takes a long time - and rowttests is a good alternative
## eg: rt = rowttests(exprs(Dilution), Dilution$liver)
\dontrun{
  data(Dilution)
  all.ttest <- ppsetApply(Dilution, ppset.ttest, covariate="liver")
}
}
\keyword{ manip }
