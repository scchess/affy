\name{expressoWidget}
\alias{expressoWidget}
\title{A widget for users to pick correction methods}
\description{
  This widget is called by expresso to allow users to select
  correction methods that will be used to process affy data.
}
\usage{
expressoWidget(BGMethods, normMethods, PMMethods, expMethods, BGDefault,
normDefault, PMDefault, expDefault) 
}
\arguments{
  \item{BGMethods}{a vector of character strings for the available methods
    that can be used as a background correction method of affy data.}
  \item{normMethods}{a vector of character strings for the available methods
    that can be used as a normalization method of affy data.}
  \item{PMMethods}{a vector of character strings for the available methods
    that can be used as a PM correction method of affy data.}
  \item{expMethods}{a vector of character strings for the available methods
    that can be used as a summary method of affy data.}
  \item{BGDefault}{a character string for the name of a default background
    correction method.}
  \item{normDefault}{a character string for the name of a default
    normalization method.}
  \item{PMDefault}{a character string for the name of a default PM correction
    method.}
  \item{expDefault}{a character string for the name of a default summary
    method.}
}
\details{
  The widget will be invoked when expresso is called with argument
  "widget" set to TRUE. Default values can be changed using the drop
  down list boxes. Double clicking on an option from the drop-down list
  makes an selection. The first element of the list for available methods
  will be the default method if no default is provided.
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

\seealso{\code{\link{expresso}}}
\examples{
if(interactive()){
  require(widgetTools)
  expressoWidget(c("mas", "none", "rma"), c("constant", "quantiles"),
c("mas", "pmonly"), c("liwong", "playerout"))
}
}
\keyword{interface}
