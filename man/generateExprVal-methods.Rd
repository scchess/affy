\name{generateExprVal}
\alias{generateExprVal}
\alias{generateExprVal-methods}
\title{Generate an expression value from the probes informations}
\description{
  Generate an expression from the from the probes
}
\usage{
generateExprVal.PPSet.container(pps.container, method, ...)
}
\arguments{
  \item{pps.container}{a \code{\link[Biobase]{container-class}} (\code{PPSet.container} for
    now.. .a \code{Plob} in the future ?..)}
  \item{method}{the method used to generate the expression value.}
  \item{...}{other parameters the method might need... (see the
    corresponding methods below...)}
}
\value{
Returns a vector of expression values.
}
\seealso{
  \code{\link{generateExprVal.method.playerout}},
  \code{\link{generateExprVal.method.avgdiff}}
  }
\examples{
  ## loading data
  data(listcel)
  data(CDF.example)

  affyid <- "A28102_at"
  pps <- get.PPSet(affyid, CDF.example, listcel[[1]])
  nprobes <- nrow(probes(pps))

  pc <- new("PPSet.container",
            pmProbes=matrix(NA, nprobes, 2, dimnames=list(NULL, rep(affyid,2))),
            mmProbes=matrix(NA, nprobes, 2)
           )

  pc[[1]] <- get.PPSet(affyid, CDF.example, listcel[[1]])
  pc[[2]] <- get.PPSet(affyid, CDF.example, listcel[[2]])

  generateExprVal(pc, method="playerout", bg.correct="bg.correct.pmonly")

}
\keyword{}
