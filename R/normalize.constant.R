normalize.AffyBatch.constant <- function(afbatch, refindex=1, FUN=mean, na.rm=TRUE) {
  
  n <- length( afbatch )
  
  if (! (refindex %in% 1:n)) stop("invalid reference index for normalization")
  refconstant <- FUN(intensity(afbatch[[refindex]]), na.rm=na.rm)
  
  #set.na.spotsd(container)
                             
  for (i in (1:n)[-refindex]) {
    m <- normalize.constant(intensity(afbatch[[i]]), refconstant, FUN=FUN, na.rm=na.rm)
    myhistory <- list(name="normalized by constant",
                      constant=attr(m,"constant"))
    attr(m,"constant") <- NULL
    intensity(afbatch)[, , i] <- m
    history(afbatch)[[i]] <- myhistory
  }
  return(afbatch)
}       



normalize.Cel.container.constant <- function(container, refindex=1, FUN=mean, na.rm=TRUE) {
  if (! inherits(container, c("Cel.container", "Cel.container.hdf5")))
    stop("container must be a 'Cel.container'")
  
  n <- length( container )
  
  if (! (refindex %in% 1:n)) stop("invalid reference index for normalization")
  refconstant <- FUN(intensity(container[[refindex]]), na.rm=na.rm)
  
  set.na.spotsd(container)
                             
  for (i in (1:n)[-refindex]) {
    m <- normalize.constant(intensity(container[[i]]), refconstant, FUN=FUN, na.rm=na.rm)
    myhistory <- list(name="normalized by constant",
                      constant=attr(m,"constant"))
    attr(m,"constant") <- NULL
    intensity(container)[, , i] <- m
    history(container)[[i]] <- myhistory
  }
  return(container)
}       



normalize.Plob.constant <- function(plob, refindex=1, FUN=mean, na.rm=TRUE) {
  ## extract the vector that will be used as a reference
  ## and apply the function FUN to generate a value
  refconstant <- FUN(c(plob@pm[, refindex], plob@mm[, refindex]), na.rm=na.rm)

  ## loop over the arrays (excluding the reference)
  for (i in (1:(nchips(plob)))[-refindex]) {
    ## only using the PM -- this is an example
    
    plob@pm[,i] <- normalize.constant(plob@pm[,i], refconstant, FUN=FUN, na.rm=na.rm)
    plob@mm[,i] <- normalize.constant(plob@mm[,i], refconstant, FUN=FUN, na.rm=na.rm)
    }
  ## state somewhere in a slot of the Plob that it has been normalized ?
  return(plob)
}


normalize.constant <- function(x, refconstant, FUN=mean, na.rm=TRUE) {
  thisconstant <- FUN(x, na.rm=na.rm)
  r <- x / thisconstant * refconstant
  attr(r,"constant") <- thisconstant * refconstant
  return(r)
}




