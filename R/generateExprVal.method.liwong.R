generateExprVal.method.liwong <- function(probes, ...) {
  probes <- t(probes)
  if (nrow(probes) == 1) {
    warning("method liwong unsuitable when only one probe pair")
    probes
  } else {
    tmp <- fit.li.wong(probes, ...)
    list(exprs=tmp$theta,se.exprs=tmp$sigma.theta)
  }
}
