## -------------------------------------------
## routine tests for expression values methods
## -------------------------------------------

library(affy)

data(affybatch.example)

i <- match("playerout", express.summary.stat.methods)
meths <- express.summary.stat.methods[-i]

for (m in meths) {
  for (mbc in pmcorrect.methods) {
     cat("expression value with method=", m, "bg correct=", mbc, "...")
     computeExprSet(affybatch.example, pmcorrect.method=mbc, summary.method=m)
     cat("done.\n")
   }
}

## playerout alone 'cause very slow
m <- "playerout"
for (mbc in pmcorrect.methods) {
  cat("expression value with method=", m, "bg correct=", mbc, "...")
  computeExprSet(affybatch.example, pmcorrect.method=mbc, summary.method=m, ids=geneNames(affybatch.example)[1:3])
  cat("done.\n")
}


