mas5 <- function(object,normalize=TRUE,sc = 500, analysis = "absolute",...){
  res <- expresso(object,bgcorrect.method="mas",pmcorrect.method="mas",normalize=FALSE,summary.method="mas",...) 
  if(normalize) 
    res <- affy.scalevalue.exprSet(res,sc=sc,analysis=analysis)
  return(res)
}

