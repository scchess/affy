\name{expressoWidget}
\alias{expressoWidget}
\title{A widget for taking user inputs for expresso}
\description{
  This widget is called by expresso of affy to allow users to input
  corection methods that will be used to process affy data.
}
\usage{
expressoWidget(BGMethods, normMethods, PMMethods, expMethods, BGDefault,
normDefault, PMDefault, expDefault) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{BGMethods}{\code{BGMethods} a vector of character strings for
    the available methods that can be used as background correction
    method of affy data}
  \item{normMethods}{\code{normMethods} a vector of character strings for
    the available methods that can be used as normalization method of
    affy data}
  \item{PMMethods}{\code{PMMethods} a vector of character strings for
    the available methods that can be used as PM correction method of
    affy data}
  \item{expMethods}{\code{expMethods}a vector of character strings for
    the available methods that can be used as summary method of affy data}
  \item{BGDefault}{\code{BGDefault} a character string for the default
    background correction method}
  \item{normDefault}{\code{normDefault} a character string for the
    default normalization method}
  \item{PMDefault}{\code{PMDefault} a character string for the default
    PM correction method}
  \item{expDefault}{\code{expDefault} a character string for the default
    summary method}
}
\details{
  Default values can be changed using the change buttons. If available
  methods are not desirable, new method may be entered using an entry
  box. However, the new method must be a valid method as the system
  currently does not check for the validity of the method entered
  through the entry box.

  The first element will be the default method if no default method is
  provided.
}
\value{
  The widget returns a list of selected correction methods.
  \item{BG}{background correction method}
  \item{NORM}{normalization method}
  \item{PM}{PM correction method}
  \item{EXP}{summary method}
}
\references{Documentations of affy package}
\author{Jianhua Zhang}
\note{The functions are part of the Bioconductor project at Dana-Farber
  Cancer Institute to provide Bioinformatic functionalities through R}

\seealso{\code{\link{tkphenoData}}, \code{link{tkSampleNames}}}
\examples{
if(interactive()){
  expressoWidget(c("111", "222", "333"), c("444", "555"), c("666",
"777"), c("888", "999"))
}
}
\keyword{misc}
