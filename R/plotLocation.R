plotLocation <- function(x, cdf=NULL, nrow=NULL, ncol=NULL, col="green", pch=22, ...) {
  if (inherits(cdf, "Cdf")) {
    nrow <- nrow(cdf@name)
    ncol <- ncol(cdf@name)
  }
  if (is.null(nrow) | is.null(ncol)) {
    stop("nrow and ncol must be specified, explicitly or with a cdf object")
  }
  points(x[,1], x[,2]
         , pch=pch, col=col, ...)
}
