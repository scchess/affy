## routine tests for the normalization methods
library(affy)

data(affybatch.example)

n.meth <- normalize.methods(affybatch.example)

## remove qspline
n.meth <- n.meth[ ! (n.meth %in% c("qspline"))]

for (m in n.meth) {
  cat("-->method=", m, "...")
  affybatch.example.n <- normalize(affybatch.example, method=m)
  cat("done.\n")
}

