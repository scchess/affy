plot.ProbeSet <- function(x, col=1:6, lty=1:6, xlab=NA, type="l", ylim=NULL, ..., covariate=NULL, col.covariate=TRUE, lty.covariate=FALSE) {
  if (is.null(ylim))
    ylim = range(c(pm(x), mm(x)), na.rm=TRUE)
  
  if (is.na(xlab))
    xlab="probes"
  
  if (! is.null(covariate)) {
    if (col.covariate)
      col <- as.integer(covariate)
    if (lty.covariate)
      lty <- as.integer(covariate)
    matplot(pm(x), col=col, lty=lty, xlab=xlab, type=type, ylim=ylim, ...)
    matplot(mm(x), col=col, lty=lty, xlab=xlab, type=type, ylim=ylim, ...)
  } else {
    matplot(pm(x), col=col, lty=lty, xlab=xlab, type=type, ylim=ylim, ...)
  }
  
}
