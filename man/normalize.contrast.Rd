\name{normalize.contrasts}
\alias{normalize.contrasts}
\alias{normalize.Cel.container.contrasts}
\alias{normalize.Plob.contrasts}
\title{Scale Cel objects in a container}
\description{
  Scale chip objects in a container.
}
\usage{
%normalize.contrast() ## currenlt maffy.normalize
normalize.Cel.container.contrasts(listcel,span=2/3,choose.subset=TRUE,subset.size=5000,verbose=TRUE,family="symmetric")
normalize.Plob.contrasts(object,span=2/3,choose.subset=TRUE,subset.size=5000,verbose=TRUE,family="symmetric")
}
\arguments{
  \item{listcel}{a \code{\link{Cel.container}}. }
  \item{object}{ a \code{\link{Plob-class}}}
  \item{span}{parameter to be passed to the function \code{\link[modreg]{loess}}.}
  \item{choose.subset}{}
  \item{subset.size}{}
  \item{verbose}{verbosity flag}
  \item{family}{parameter to be passed to the function \code{\link[modreg]{loess}}.}
}
\value{
  An object of the same class as the one passed.
}
\author{
  L. Gautier <laurent@cbs.dtu.dk> (<--- not me... this is a copy/paste
  accident. I think the author is Rafael or Ben Bolstad...)
}
\examples{
    data(listcel)
    data(CDF.example)
    plob <- convert(listcel,CDF.example)
    tmp <- normalize(listcel,CDF.example,method="contrasts",subset.size=1000)
    tmp <- normalize(plob,method="contrasts",subset.size=1000)
}

\seealso{
  \code{\link{maffy.normalize}}
  }
%\examples{
%  
%}
\keyword{manip}
