normalize.AffyBatch.contrasts <- function(abatch,span=2/3,choose.subset=TRUE,subset.size=5000,verbose=TRUE,family="symmetric",pmonly=FALSE) {
  
  if(pmonly)
    Index <- unlist(pmindex(abatch))
  else
    Index <- c(unlist(pmindex(abatch)),unlist(mmindex(abatch)))
  
  
  ##we need default argumetns becuase they are used in this transitional file
  alldata <- intensity(abatch)[Index,]
  
  if(choose.subset)
    subset1 <- maffy.subset(alldata,verbose=verbose,subset.size=subset.size)$subset
  else
    subset1 <- sample(1:dim(alldata)[1],subset.size)
  aux <-   maffy.normalise(alldata,subset=subset1,verbose=verbose,span=span,family=family)
  
  intensity(abatch)[Index,] <- aux
  
  return(abatch)
}

normalize.Cel.container.contrasts <- function(listcel,span=2/3,choose.subset=TRUE,subset.size=5000,verbose=TRUE,family="symmetric") { 
  
  cols <- length(listcel)
  rows <- length(intensity(listcel[[1]]))
  chipdim <- dim(intensity(listcel[[1]]))
  
  x <- matrix(0,rows,cols)
  for (i in 1:cols) x[,i] <- c(intensity(listcel[[i]]))
  
  n <-  dim(x)[1]/2
  if(choose.subset)
    subset1 <- maffy.subset(x,verbose=verbose,subset.size=subset.size)$subset
  else
    subset1 <- sample(1:dim(x)[1],subset.size)

  x <-   maffy.normalise(x,subset=subset1,verbose=verbose,span=span,family=family)

  set.na.spotsd(listcel) # set 'sd' to nothing (meaningless after normalization)
  
  for (i in 1:cols) {
    intensity(listcel)[, , i]  <- matrix(x[,i],chipdim[1], chipdim[2])
    history(listcel)[[i]]$name <- "normalized by contrasts"
  }
  return(listcel)
}








