generateExprVal.method.liwong <- function(probes, ...) {
  matos <- as.data.frame(lapply(probes, function(x) {x$pm - x$mm}))
  matos <- t(as.matrix(matos))
  fit.li.wong(matos, ...)$theta
}
