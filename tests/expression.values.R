## routine tests for expression values methods
library(affy)


## Cel.container and Cdf
cat("Cel.container and Cdf:\n")

data(CDF.example)
data(listcel)

## "playerout" very slow. tested individually.
i <- match("playerout", express.summary.stat.methods)
express.summary.stat.methods <- express.summary.stat.methods[-i]

for (m in express.summary.stat.methods) {
  for (mbc in bg.correct.methods) {
    cat("expression value with method=", m, "bg correct=", mbc, "...")
    generateExprSet(listcel, CDF.example, method=m, bg.correct=mbc)
    cat("done.\n")
  }
}

m <- "playerout"
for (mbc in bg.correct.methods) {
  cat("expression value with method=", m, "bg correct=", mbc, "...")
  generateExprSet(listcel, CDF.example, method=m, bg.correct=mbc, ids=CDF.example@name.levels[1:3])
  cat("done.\n")
}

## Plob

cat("Plob:\n")
data(Plob)

#express()
