generateExprVal.method.liwong <- function(probes, ...) {
  probes <- t(probes)
  if (nrow(probes) == 1) {
    warning("method liwong unsuitable when only one probe pair")
    probes
  } else {
    fit.li.wong(probes, ...)$theta
  }
}
