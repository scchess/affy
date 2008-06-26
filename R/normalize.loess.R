normalize.AffyBatch.loess <- function(abatch,type=c("together","pmonly","mmonly","separate"),...) {

  type <- match.arg(type)

  if (type == "separate"){
    Index <- unlist(indexProbes(abatch,"pm"))
    intensity(abatch)[Index,] <- normalize.loess(intensity(abatch)[Index,], ...)
    Index <- unlist(indexProbes(abatch,"mm"))
    intensity(abatch)[Index,] <- normalize.loess(intensity(abatch)[Index,], ...)
  } else if (type=="together"){
    Index <- unlist(indexProbes(abatch,"both"))
    intensity(abatch)[Index,] <- normalize.loess(intensity(abatch)[Index,], ...)
  } else if (type=="pmonly"){
    Index <- unlist(indexProbes(abatch,"pm"))
    intensity(abatch)[Index,] <- normalize.loess(intensity(abatch)[Index,], ...)
  } else if (type=="mmonly"){
    Index <- unlist(indexProbes(abatch,"mm"))
    intensity(abatch)[Index,] <- normalize.loess(intensity(abatch)[Index,], ...)
  }
  ##set.na.spotsd(listcel) # set 'sd' to nothing (meaningless after normalization)
  ##cat(cols,rows)


  ##need to use MIAME
  ##for (i in 1:abatch@nexp) {
  ##  history(abatch)[[i]] <- list(name="normalized by loess")
  ##}

  return(abatch)
}



normalize.loess <- function(mat, subset=sample(1:(dim(mat)[1]), min(c(5000, nrow(mat)))),
                            epsilon=10^-2, maxit=1, log.it=TRUE, verbose=TRUE, span=2/3,
                            family.loess="symmetric"){

  J <- dim(mat)[2]
  II <- dim(mat)[1]
  if(log.it){
    mat <- log2(mat)
  }

  change <- epsilon +1
  iter <- 0
  w <- c(0, rep(1,length(subset)), 0) ##this way we give 0 weight to the
                                      ##extremes added so that we can interpolate

  while(iter < maxit){
    iter <- iter + 1
    means <- matrix(0,II,J) ##contains temp of what we substract

    for (j in 1:(J-1)){
      for (k in (j+1):J){
        y <- mat[,j] - mat[,k]
        x <- (mat[,j] + mat[,k]) / 2
        index <- c(order(x)[1], subset, order(-x)[1])
        ##put endpoints in so we can interpolate
        xx <- x[index]
        yy <- y[index]
        aux <-loess(yy~xx, span=span, degree=1, weights=w, family=family.loess)
        aux <- predict(aux, data.frame(xx=x)) / J
        means[, j] <- means[, j] + aux
        means[, k] <- means[, k] - aux
        if (verbose)
          cat("Done with",j,"vs",k,"in iteration",iter,"\n")
      }
    }
    mat <- mat - means
    change <- max(colMeans((means[subset,])^2))

    if(verbose)
      cat(iter, change,"\n")

  }

  if ((change > epsilon) & (maxit > 1))
    warning(paste("No convergence after", maxit, "iterations.\n"))

  if(log.it) {
    return(2^mat)
  } else
    return(mat)
}
