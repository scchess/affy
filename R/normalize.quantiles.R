normalize.AffyBatch.quantiles <- function(abatch,pmonly=FALSE) {

  pms <- unlist(pmindex(abatch))
  noNA <- apply(intensity(abatch)[pms,],1,function(x) all(!is.na(x)))
  pms <- pms[noNA]
  intensity(abatch)[pms,] <- normalize.quantiles(intensity(abatch)[pms, ])
  if(!pmonly){ 
    mms <- unlist(mmindex(abatch))
    noNA <- apply(intensity(abatch)[mms,],1,function(x) all(!is.na(x)))
    mms <- mms[noNA]

    intensity(abatch)[mms,] <- normalize.quantiles(intensity(abatch)[mms, ])
  }
  
  ##this is MIAME we need to decide how to do this properly.
  ##for (i in 1:length(abatch)) {
  ##  history(abatch)[[i]]$name <- "normalized by quantiles"
  ##}

                return(abatch)
}
  
normalize.quantiles <- function(x){

  rows <- dim(x)[1]
  cols <- dim(x)[2]

  if (!is.matrix(x)){
    stop("Matrix expected in normalize.quantiles")
  }
  
  matrix(.C("qnorm_c", as.double(as.vector(x)), as.integer(rows), as.integer(cols))[[1]], rows, cols)
}


normalize.AffyBatch.quantiles.robust <- function(abatch, pmonly=FALSE) {

  pms <- unlist(pmindex(abatch))
  intensity(abatch)[pms, ] <- normalize.quantiles.robust(intensity(abatch)[pms, ])
  if(!pmonly){ 
    mms <- unlist(mmindex(abatch))
    intensity(abatch)[mms, ] <- normalize.quantiles.robust(intensity(abatch)[mms, ])
  }
  
  ##this is MIAME we need to decide how to do this properly.
  ##for (i in 1:length(abatch)) {
  ##  history(abatch)[[i]]$name <- "normalized by quantiles"
  ##}

  return(abatch)
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
