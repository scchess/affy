## routine tests for expression values methods
library(affy)


## Cel.container and Cdf
cat("Cel.container and Cdf:\n")

data(CDF.example)
data(listcel)

for (m in express.summary.stat.methods) {
  for (mbc in bg.correct.methods) {
    cat("expression value with method=", m, "bg correct=", mbc, "...")
    generateExprSet(listcel, CDF.example, method=m, bg.correct=mbc)
    cat("done.\n")
  }
}

## Plob

cat("Plob:\n")
data(Plob)

#express()
