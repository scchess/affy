normalize.Plob.contrasts <- function(object,span=2/3,choose.subset=T,subset.size=5000,verbose=T,family="symmetric") {
  ##we need default argumetns becuase they are used in this transitional file
  alldata <- rbind(object@pm,object@mm)
  n <-  dim(alldata)[1]/2
  if(choose.subset)
    subset1 <- maffy.subset(alldata,verbose=verbose,subset.size=subset.size)$subset
  else
    subset1 <- sample(1:dim(alldata)[1],subset.size)
  aux <-   maffy.normalise(alldata,subset=subset1,verbose=verbose,span=span,family=family)

  
  pm(object) <- aux[1:n,]
  mm(object) <- aux[(n+1):(2*n),]

  
  return(object)
}

normalize.Cel.contrasts <- function(listcel,span=2/3,choose.subset=T,subset.size=5000,verbose=T,family="symmetric") { 
  
  cols <- length(listcel)
  rows <- length(listcel[[1]]@intensity)
  chipdim <- dim(listcel[[1]]@intensity)
  
  x <- matrix(0,rows,cols)
  for (i in 1:cols) x[,i] <- c(listcel[[i]]@intensity)
  
  n <-  dim(x)[1]/2
  if(choose.subset)
    subset1 <- maffy.subset(x,verbose=verbose,subset.size=subset.size)$subset
  else
    subset1 <- sample(1:dim(x)[1],subset.size)

  x <-   maffy.normalise(x,subset=subset1,verbose=verbose,span=span,family=family)
  
  for (i in 1:cols) {
    listcel[[i]]@intensity  <- matrix(x[,i],chipdim[1], chipdim[2])
    listcel[[i]]@history$name <- "normalized by quantiles"
    listcel[[i]]@sd <- matrix() # set 'sd' to nothing (meaningless after normalization)
  }
  return(listcel)
}
