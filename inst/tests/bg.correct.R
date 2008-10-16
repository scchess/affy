library(affy)

data(affybatch.example)

meth <- bgcorrect.methods

cat("background correction:\n")

for (m in meth) {
  cat(m,"...")
  abatch.bgc <- bg.correct(affybatch.example, method=m)
  cat("done.\n")
}
