## Currently, the input is a 2 matrices a pm and a mm

generateExprVal.method.avgdiff <- function(probes, ...) {
  apply(probes, 2, mean)
}
