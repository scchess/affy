generateExprVal.method.liwong <- function(probes, ...) {
  matos <- as.data.frame(lapply(probes, function(x) {x$pm - x$mm}))
  matos <- t(as.matrix(matos))
  if (nrows(matos) == 1) {
    warning("method liwong unsuitable when only one probe pair")
    matos
  } else {
    fit.li.wong(matos, ...)$theta
  }
}
