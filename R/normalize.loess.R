normalize.Cel.container.loess <- function(listcel, ...) {
  
  cols <- length(listcel)
  rows <- length(intensity(listcel[[1]])) # assuming all the Cel are of the same size in listcel
  chipdim <- dim(intensity(listcel[[1]]))
  
  ## this is may be too much.. I am taking every intensity in there... something like
  ## taking only the PM (see invariantset) could be what you want (this is only a trial,
  ## know better your own algorithm that I will ever... modify if needed).  
  x <- matrix(0,rows,cols)
  for (i in 1:cols) x[,i] <- c(intensity(listcel[[i]]))

  x <- normalize.loess(x, ...)

  set.na.spotsd(listcel) # set 'sd' to nothing (meaningless after normalization)
  ##cat(cols,rows)
  for (i in 1:cols) {
    intensity(listcel)[, , i]  <- matrix(x[,i],chipdim[1], chipdim[2])
    history(listcel)[[i]] <- list(name="normalized by loess")
  }

  return(listcel)
}

normalize.Plob.loess <- function(plob, ...){

  x <-  matrix(0, nprobes(plob)*2, nchips(plob))
  
  for (i in 1:nchips(plob)) {
    x[1:nprobes(plob), i] <- pm(plob)[, i]
    x[(nprobes(plob)+1):(nprobes(plob)*2), i] <- mm(plob)[, i]
  }

  y <- normalize.loess(x, ...)
  
  for (i in 1:nchips(plob)) {
    pm(plob)[, i] <- y[1:nprobes(plob), i]
    mm(plob)[, i] <- y[(nprobes(plob)+1):(nprobes(plob)*2), i]
  }
  
  return(plob)
  
}


normalize.loess <- function(mat, subset=sample(1:(dim(mat)[1]), min(c(5000, nrow(mat)))),
                            epsilon=10^-2, maxit=1, log.it=TRUE, verbose=TRUE, span=2/3,
                            family.loess="symmetric"){
  
  J <- dim(mat)[2]
  II <- dim(mat)[1]
  newData <- mat
  if(log.it){
    mat <- log2(mat)
    newData <- log2(newData)
  }
  
  change <- epsilon +1
  fs <- matrix(0, II, J)##contains what we substract
  iter <- 0
  w <- c(0, rep(1,length(subset)), 0) ##this way we give 0 weight to the
                                      ##extremes added so that we can interpolate
  
  while(iter < maxit){
    iter <- iter + 1
    means <- matrix(0,II,J) ##contains temp of what we substract
    
    for (j in 1:(J-1)){
      for (k in (j+1):J){
        y <- newData[,j] - newData[,k]
        x <- (newData[,j] + newData[,k]) / 2
        index <- c(order(x)[1], subset, order(-x)[1])
        ##put endpoints in so we can interpolate
        xx <- x[index]
        yy <- y[index]
        aux <-loess(yy~xx, span=span, degree=1, weights=w, family=family.loess)
        aux <- predict.loess(aux, data.frame(xx=x)) / J
        means[, j] <- means[, j] + aux 
        means[, k] <- means[, k] - aux
        if (verbose)
          cat("Done with",j,"vs",k," in iteration ",iter,"\n")
      }
    }
    fs <- fs + means
    newData <- mat - fs
    change <- max(apply((means[subset,])^2, 2, mean))
    
    if(verbose)
      cat(iter, change,"\n")
    
    oldfs <- fs
    
  }
  
  if ((change > epsilon) & (maxit > 1))
    warning(paste("No convergence after", maxit, "iterations.\n"))
  
  if(log.it) {
    return(2^newData)
  } else
    return(newData)
}
