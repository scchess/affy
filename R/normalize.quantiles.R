normalize.AffyBatch.quantiles <- function(abatch,pmonly=FALSE) {

  pms <- unlist(pmindex(abatch))
  intensity(abatch)[pms,] <- normalize.quantiles(intensity(abatch)[pms,])
  if(!pmonly){ 
    mms <- unlist(mmindex(abatch))
    intensity(abatch)[mms,] <- normalize.quantiles(intensity(abatch)[mms,])
  }
  
  ##this is MIAME we need to decide how to do this properly.
  ##for (i in 1:length(abatch)) {
  ##  history(abatch)[[i]]$name <- "normalized by quantiles"
  ##}

  return(abatch)
}

normalize.Cel.container.quantiles <- function(listcel) {
  ## make the matrix (not very memory friendly... will be better when HDF5,
  ## or an another mecanism is around...)
  ## note to Rafael: I hope I did not butcher your code too much... this is just
  ## an example on how to get the function to the new objects scheme
  ## did with 'Cel' and the corresponding generic method (see 'Cel.R'), and of course
  ##
  ## Laurent
  
  rows <- length(intensity(listcel)[, , 1])
  cols <- length(listcel) # assuming all the Cel are of the same size in listcel
  chipdim <- dim(intensity(listcel)[, , 1])
  
  ## this is may be too much.. I am taking every intensity in there... something like
  ## taking only the PM (see invariantset) could be what you want (this is only a trial,
  ## know better your own algorithm that I will ever... modify if needed).
  
  x <- matrix(0, rows, cols)
  
  for (i in 1:cols)
    x[, i] <- c(intensity(listcel)[, , i])

  x <- normalize.quantiles(x)
  set.na.spotsd(listcel) # set 'sd' to nothing (meaningless after normalization)
  for (i in 1:cols) {
    intensity(listcel)[, , i]  <- matrix(x[,i],chipdim[1], chipdim[2])
    history(listcel)[[i]]$name <- "normalized by quantiles"
  }

  return(listcel)
}


normalize.quantiles <- function(x){

  rows <- dim(x)[1]
  cols <- dim(x)[2]
  
  matrix(.C("qnorm_c", as.double(as.vector(x)), as.integer(rows), as.integer(cols))[[1]], rows, cols)
}


normalize.quantiles.robust <- function(x,weights=NULL,remove.extreme=c("variance","mean","both","none"),n.remove=1,approx.meth = FALSE,...){
  
  calc.var.ratios <- function(x){
    cols <- dim(x)[2]
    vars <- apply(x,2,var)
    results <- matrix(0,cols,cols)
    for (i in 1:cols-1)
      for (j in (i+1):cols){
        results[i,j] <- vars[i]/vars[j]
        results[j,i] <- vars[j]/vars[i]
      }
    results
  }

  calc.mean.dists <- function(x){
    cols <- dim(x)[2]
    means <- apply(x,2,mean)
    results <- matrix(0,cols,cols)
    for (i in 1:cols-1)
      for (j in (i+1):cols){
        results[i,j] <- means[i] - means[j]
        results[j,i] <- means[j] - means[i]
      }
    results
  }
  
  rows <- dim(x)[1]
  cols <- dim(x)[2]
  
  if (is.null(weights)){
    weights <- rep(1,cols)
    if (remove.extreme == "variance"){
      var.ratios <- calc.var.ratios(x)
      vars.big <- apply(var.ratios,1,sum)
      vars.small <- apply(var.ratios,2,sum)
      var.adj <- vars.big + vars.small
      remove.order <- order(-var.adj)
      weights[remove.order[1:n.remove]] <- 0
    }
    if (remove.extreme == "mean"){
      means <- abs(apply(calc.mean.dists(x),2,sum))
      remove.order <- order(-means)
      weights[remove.order[1:n.remove]] <- 0
    }
    if (remove.extreme == "both"){
      var.ratios <- calc.var.ratios(x)
      vars.big <- apply(var.ratios,1,sum)
      vars.small <- apply(var.ratios,2,sum)
      var.adj <- vars.big + vars.small
      means <- abs(apply(calc.mean.dists(x),2,sum))
      # by convention we will remove first the most extreme variance, then the most extreme mean
      remove.order <- order(-var.adj)
      weights[remove.order[1]] <- 0
      remove.order <- order(-means)
      weights[remove.order[1]] <- 0
    }
  }
  if (length(weights) != cols){
    stop("Weights vector incorrect length\n")
  }
  if (sum(weights > 0) < 2){
    stop("Need at least two non negative weights\n")
  }
  cat("Chip weights are ",weights,"\n") 
  if (approx.meth == FALSE){
    matrix(.C("qnorm_robust_c",as.double(as.vector(x)),as.double(weights),as.integer(rows),as.integer(cols))[[1]],rows,cols)
  } else {
    cat("Approximation currently not implemented \nFalling back to standard Quantile method\n")
    matrix(.C("qnorm_robust_c",as.double(as.vector(x)),as.double(weights),as.integer(rows),as.integer(cols))[[1]],rows,cols)
  }
}
