normalize.AffyBatch.invariantset <- function(abatch, prd.td=c(0.003,0.007), progress=FALSE) {

  require(modreg, quietly=TRUE)
  
  w.pm <- unlist(indexProbes(abatch, which="pm"))             # boolean to find the PM probes
  i.pm <- rep(FALSE, nrow(abatch) * ncol(abatch))
  i.pm[w.pm] <- TRUE
  rm(w.pm)
  
  np <- sum(i.pm)                                     # number of PM probes
  nc  <-  length(abatch)                                 # number of CEL files
  
  # take as a reference the array having the median overall intensity
  m <- vector("numeric", length=nc)
  for (i in 1:nc)
    m[i] <- mean(intensity(abatch)[, i][i.pm])
  refindex <- trunc(median(rank(m)))
  rm(m)           

  if (progress) cat("Data from", chipNames(abatch)[refindex], "used as baseline.\n")
  
  ##set.na.spotsd(cel.container)
  
  normhisto <- vector("list", length=nc)
  normhisto[[refindex]] <- list(name="reference for the invariant set")
  
  ## loop over the CEL files and normalize them
  for (i in (1:nc)[-refindex]) {
  
    if (progress) cat("normalizing array", chipNames(abatch)[i], "...")
    
    ##temporary
    tmp <- normalize.invariantset(c(intensity(abatch)[, i])[i.pm],
                                  c(intensity(abatch)[, refindex])[i.pm],
                                  prd.td)
    i.set <- which(i.pm)[tmp$i.set]
    tmp <- as.numeric(approx(tmp$n.curve$y, tmp$n.curve$x,
                             xout=intensity(abatch)[, i], rule=2)$y)
    attr(tmp,"invariant.set") <- NULL
    intensity(abatch)[, i] <- tmp

    ## storing information about what has been done
    normhisto[[i]] <- list(name="normalized by invariant set",
                           invariantset=i.set)
    
    if (progress) cat("done.\n")
    
  }
  
  attr(abatch, "normalization") <- normhisto
  return(abatch)
}



##  The 'common-to-all' part of the algorithm. Operates on two vectors of numeric data
##
normalize.invariantset <- function(data, ref, prd.td=c(0.003,0.007)) {

  np <- length(data)
  r.ref <- rank(ref)
  r.array <- rank(data)
  
  ## init
  prd.td.adj <- prd.td*10                           # adjusted threshold things
  i.set <- rep(TRUE, np)                            # index all the PM probes as being in the invariant set
  ns <- sum(i.set)                                  # number of probes in the invariant set
  ns.old <- ns+50+1                                 # number of probes previously in the invariant set
    
  ## iterate while the number of genes in the invariant set (ns) still varies...
  while ( (ns.old-ns) > 50 ) {
    air <- (r.ref[i.set] + r.array[i.set]) / (2*ns)  # average intensity rank for the probe intensities
    prd <- abs(r.ref[i.set] - r.array[i.set]) / ns
    threshold <- (prd.td.adj[2]-prd.td[1]) * air + prd.td.adj[1]
    i.set[i.set] <- (prd < threshold)
    
    ns.old <- ns
    ns <- sum(i.set)
    
    if (prd.td.adj[1] > prd.td[1])
      prd.td.adj <- prd.td.adj * 0.9  # update the adjusted threshold parameters
  }
  
  ## the index i.set corresponds to the 'invariant genes'
  n.curve <- smooth.spline(ref[i.set], data[i.set])
  ## n.curve$x contains smoothed reference intensities
  ## n.curve$y contains smoothed i-th array intensities
  
  ##data <- as.numeric(approx(n.curve$y, n.curve$x, xout=data)$y)
  ##attr(data,"invariant.set") <- i.set
  ##return(data)
  return(list(n.curve=n.curve, i.set=i.set))
}







