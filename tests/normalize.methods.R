## routine tests for the normalization methods
library(affy)

## Cel.container and Cdf
cat("Cel.container and Cdf\n")

data(CDF.example)
data(listcel)
for (m in normalize.methods(listcel)) {
  cat("-->", m, "...")
  listcel.n <- normalize(listcel, CDF.example, method=m)
  cat("done.\n")
}

## Plob
cat("Plob\n")

data(Dilution)
for (m in normalize.methods(listcel)) {
  cat("-->", m, "...")
  listcel.n <- normalize(Dilution, method=m)
  cat("done.\n")
}
