generateExprVal.method.mas <- function(probes, ...)
{
  
  probes <- log2(probes)
  
  slg <- rep(NA, ncol(probes))
  
  for (i in 1:ncol(probes)) {
    
    slg[i] <- tukey.biweight(probes[ ,i], ...)
    
  }
  
  return(slg)
    
}

affy.scalevalue.exprSet <- function(eset, sc=500, analysis)
{
  
  analysis <- match(analysis, c("absolute", "comparison"))
  
  if (analysis == "absolute")
    nf <- 1
  else
    stop("not implemented")
  
  for (i in 1:ncol(exprs(eset))) {
    slg <- exprs(eset)[, i]
    sf <- sc / mean(2^slg, trim=0.02)  
    reported.value <- nf * sf * 2^slg
    eset@exprs[, i] <- reported.value
  }
  
  return(eset)
}
