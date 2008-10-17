## routine tests for the normalization methods
library(affy)

data(Dilution)

n.meth <- normalize.methods(Dilution)

## remove qspline
##n.meth <- n.meth[ ! (n.meth %in% c("qspline"))]

for (m in n.meth) {
  cat("-->method=", m, "...")
  affybatch.example.n <- normalize(Dilution, method=m)
  cat("done.\n")
}
