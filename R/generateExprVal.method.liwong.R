generateExprVal.method.liwong <- function(matos, ...) {
  matos <- t(matos)
  if (nrow(matos) == 1) {
    warning("method liwong unsuitable when only one probe pair")
    matos
  } else {
    fit.li.wong(matos, ...)$theta
  }
}
