normalize.AffyBatch.contrasts <- function(abatch,span=2/3,choose.subset=TRUE,subset.size=5000,verbose=TRUE,family="symmetric",pmonly=FALSE) {
  
  if(pmonly)
    Index <- unlist(pmindex(abatch))
  else
    Index <- unlist(indexProbes(abatch,"both"))
  
  
  ##we need default argumetns becuase they are used in this transitional file
  alldata <- intensity(abatch)[Index,]
  
  if(choose.subset)
    subset1 <- maffy.subset(alldata,verbose=verbose,subset.size=subset.size)$subset
  else
    subset1 <- sample(1:dim(alldata)[1],subset.size)
  aux <-   maffy.normalize(alldata,subset=subset1,verbose=verbose,span=span,family=family)
  
  intensity(abatch)[Index,] <- aux

  ##attr(abatch, "normalization") <- normhisto
  return(abatch)
}









