## routine tests for expression values methods
library(affy)

## Cel.container and Cdf
cat("Cel.container and Cdf:\n")

data(CDF.example)
data(listcel)

## "playerout" very slow. tested individually.
i <- match("playerout", express.summary.stat.methods)
meths <- express.summary.stat.methods[-i]

for (m in meths) {
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


## AffyBatch

cat("AffyBatch:\n")
data(affybatch.example)

i <- match("playerout", express.summary.stat.methods)
meths <- express.summary.stat.methods[-i]

for (m in meths) {
  for (mbc in bg.correct.methods) {
     cat("expression value with method=", m, "bg correct=", mbc, "...")
     computeExprSet(affybatch.example, bg.method=mbc, summary.method=m, warnings=FALSE)
     cat("done.\n")
   }
}

m <- "playerout"
for (mbc in bg.correct.methods) {
  cat("expression value with method=", m, "bg correct=", mbc, "...")
  computeExprSet(affybatch.example, bg.method=mbc, summary.method=m, ids=geneNames(affybatch.example)[1:3])
  cat("done.\n")
}


#express()
