## (see normalize.invariantset.R to know more....)

normalize.Cel.container.constant <- function(container, refindex=1, FUN=mean, na.rm=TRUE) {
  if (! inherits(container, "Cel.container"))
    stop("container must be a 'Cel.container'")
  
  n <- length( container )
  
  if (! (refindex %in% 1:n)) stop("invalid reference index for normalization")
  
  refconstant <- FUN(container[[refindex]]@intensity, na.rm=na.rm)
  
  for (i in 1:n) {
    container[[i]]@intensity <- normalize.constant(container[[i]]@intensity, refconstant, FUN=FUN, na.rm=na.rm)
    container[[i]]@history <- list(name="normalized by constant",
                                 constant=attr(container[[i]]@intensity,"constant"))
    container[[i]]@sd <- matrix()
    attr(container[[i]]@intensity, "constant") <- NULL
  }
  return(container)
}       


normalize.Plob.constant <- function(plob, refindex=1, FUN=mean, na.rm=TRUE) {
  ## extract the vector that will be used as a reference
  ## and apply the function FUN to generate a value
  refconstant <- FUN(c(plob@pm[, refindex], plob@mm[, refindex]), na.rm=na.rm)

  ## loop over the arrays (excluding the reference)
  for (i in (1:plob@nchips)[-refindex]) {
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




