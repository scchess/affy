##################################################################
##
## file: normalize.quantiles.R
##
## For a description of quantile normalization method see
##
##  Bolstad, B. M., Irizarry R. A., Astrand, M, and Speed, T. P. (2003)(2003)
##  A Comparison of Normalization Methods for High
##  Density Oligonucleotide Array Data Based on Bias and Variance.
##  Bioinformatics 19,2,pp 185-193
##
## History
## Pre Aug 23, 2003 Two years worth of stuff
## Aug 23, 2003 - Added use.log2 to "robust",
##                added ability to pass additional parameters
##                to normalize.AffyBatch.Quantiles.robust
##                changed pmonly parameters on functions
##                so that it is now a string argument "type"
##                the options are pmonly, mmonly, together, separate
## Jan 31, 2004 - put a check for an integer matrix and force coercision to
##                doubles if required in normalize.quantiles
## Mar 13, 2005 - Modifications to normalize.quantiles.robust including removing
##                approx.method which never got implemented. Making it a use a .Call()
##                rather than a .C()
##
## Sep 20, 2006 - fix .Call in normalize.quantiles.robust
##
##################################################################

normalize.AffyBatch.quantiles <- function(abatch,type=c("separate","pmonly","mmonly","together")) {


  type <- match.arg(type)

  if ((type == "pmonly")|(type == "separate")){
    pms <- unlist(pmindex(abatch))
    ## Change to faster computation of noNA - SDR 11/06/2003
    ##noNA <- apply(intensity(abatch)[pms,,drop=FALSE],1,function(x) all(!is.na(x)))
    noNA <- rowSums(is.na(intensity(abatch)[pms,,drop=FALSE])) == 0
    pms <- pms[noNA]
    intensity(abatch)[pms,] <- normalize.quantiles(intensity(abatch)[pms,,drop=FALSE ],copy=FALSE)
  }
  if((type == "mmonly") | (type == "separate")){
    mms <- unlist(mmindex(abatch))
    ## Change to faster computation of noNA - SDR 11/06/2003
    ##noNA <- apply(intensity(abatch)[mms,,drop=FALSE],1,function(x) all(!is.na(x)))
    noNA <- rowSums(is.na(intensity(abatch)[mms,,drop=FALSE])) == 0
    mms <- mms[noNA]

    intensity(abatch)[mms,] <- normalize.quantiles(intensity(abatch)[mms,,drop=FALSE ],copy=FALSE)
  }
  if (type == "together"){
    pms <- unlist(indexProbes(abatch,"both"))
    intensity(abatch)[pms,]  <- normalize.quantiles(intensity(abatch)[pms,,drop=FALSE ],copy=FALSE)
  }

  ##this is MIAME we need to decide how to do this properly.
  ##for (i in 1:length(abatch)) {
  ##  history(abatch)[[i]]$name <- "normalized by quantiles"
  ##}

  return(abatch)
}

normalize.quantiles <- function(x,copy=TRUE){

  rows <- dim(x)[1]
  cols <- dim(x)[2]

  if (!is.matrix(x)){
    stop("Matrix expected in normalize.quantiles")
  }

  if (is.integer(x)){
    x <- matrix(as.double(x),rows,cols)
    copy <- FALSE
  }

  #matrix(.C("qnorm_c", as.double(as.vector(x)), as.integer(rows), as.integer(cols))[[1]], rows, cols)

  .Call("R_qnorm_c",x,copy, PACKAGE="affy");
}


normalize.AffyBatch.quantiles.robust <- function(abatch, type=c("separate","pmonly","mmonly","together"),weights=NULL,remove.extreme=c("variance","mean","both","none"),n.remove=1,use.median=FALSE,use.log2=FALSE) {

  type <- match.arg(type)

  if ((type == "pmonly")|(type == "separate")){
    pms <- unlist(pmindex(abatch))
    intensity(abatch)[pms, ] <- normalize.quantiles.robust(intensity(abatch)[pms, ], copy=FALSE,weights=weights,remove.extreme,n.remove=n.remove,use.median=use.median,use.log2=use.log2)
  }
  if ((type == "mmonly")|(type == "separate")){
    mms <- unlist(mmindex(abatch))
    intensity(abatch)[mms, ] <- normalize.quantiles.robust(intensity(abatch)[mms, ],copy=FALSE,weights=weights,remove.extreme,n.remove=n.remove,use.median=use.median,use.log2=use.log2)
  }

  if (type == "together"){
    pms <- unlist(indexProbes(abatch,"both"))
    intensity(abatch)  <- normalize.quantiles.robust(intensity(abatch)[pms,,drop=FALSE ],copy=FALSE, weights=weights,remove.extreme=remove.extreme,n.remove=n.remove,use.median=use.median,use.log2=use.log2)
  }



  ##this is MIAME we need to decide how to do this properly.
  ##for (i in 1:length(abatch)) {
  ##  history(abatch)[[i]]$name <- "normalized by quantiles"
  ##}

  return(abatch)
}

normalize.quantiles.robust <- function(x,copy=TRUE,weights=NULL,remove.extreme=c("variance","mean","both","none"),n.remove=1,use.median=FALSE,use.log2=FALSE){

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
    means <- colMeans(x)
    results <- matrix(0,cols,cols)
    for (i in 1:cols-1)
      for (j in (i+1):cols){
        results[i,j] <- means[i] - means[j]
        results[j,i] <- means[j] - means[i]
      }
    results
  }

  use.huber <- FALSE
  remove.extreme <- match.arg(remove.extreme)

  rows <- dim(x)[1]
  cols <- dim(x)[2]

  if (is.null(weights)){
    weights <- rep(1,cols)
    if (remove.extreme == "variance"){
      var.ratios <- calc.var.ratios(x)
      vars.big <- rowSums(var.ratios)
      vars.small <- colSums(var.ratios)
      var.adj <- vars.big + vars.small
      remove.order <- order(-var.adj)
      weights[remove.order[1:n.remove]] <- 0
    }
    if (remove.extreme == "mean"){
        means <- abs(colSums(calc.mean.dists(x)))
      remove.order <- order(-means)
      weights[remove.order[1:n.remove]] <- 0
    }
    if (remove.extreme == "both"){
      var.ratios <- calc.var.ratios(x)
      vars.big <- rowSums(var.ratios)
      vars.small <- colSums(var.ratios)
      var.adj <- vars.big + vars.small
      means <- abs(colSums(calc.mean.dists(x)))
      # by convention we will remove first the most extreme variance, then the most extreme mean
      remove.order <- order(-var.adj)
      weights[remove.order[1]] <- 0
      remove.order <- order(-means)
      weights[remove.order[1]] <- 0
    }
  } else {

    if (is.numeric(weights)){
      if (length(weights) != cols){
        stop("Weights vector incorrect length\n")
      }
      if (sum(weights > 0) < 1){
        stop("Need at least one non negative weights\n")
      }
      if (any(weights < 0)){
        stop("Can't have negative weights")
      }

      
    } else {
      if (weights =="huber"){
        use.huber <- TRUE
        weights <- rep(1,cols)
      } else {
        stop("Don't recognise weights argument as valid.")
      }

    }

  }


      
  cat("Chip weights are ",weights,"\n")
  ###matrix(.C("qnorm_robust_c",as.double(as.vector(x)),as.double(weights),as.integer(rows),as.integer(cols),as.integer(use.median),as.integer(use.log2),as.integer(use.huber),
  ####PACKAGE="affy")[[1]],rows,cols)
  
  ####R_qnorm_robust_c(SEXP X, SEXP copy, SEXP R_weights, SEXP R_use_median, SEXP R_use_log2, SEXP R_weight_scheme)
  .Call("R_qnorm_robust_c",x,copy,weights,as.integer(use.median),as.integer(use.log2),as.integer(use.huber),PACKAGE="affy")


  
}
