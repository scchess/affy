plot.ProbeSet <- function(x, which=c("pm", "mm"), xlab="probes", type="l", ylim=NULL, ...) {

  which <- match.arg(which)
  ## Both should be f. Also, use the generic directly. SDR 10/21/03
#   if (which == "pm")
#     f <- getMethod("pm", "ProbeSet")
#   else
#     g <- getMethod("mm", "ProbeSet")
  f <- if (which == "pm")
      pm
  else mm
  
  if (is.null(ylim))
    ylim = range(c(f(x)), na.rm=TRUE)
  
  if (is.na(xlab))
    xlab="probes"
  
  matplot(f(x), xlab=xlab, type=type, ylim=ylim, ...)  
}
