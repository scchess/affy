## (see normalize.invariantset.R to know more....)

normalize.Cel.constant <- function(container, refindex=1) {
  if (! inherits(container, "Cel.container"))
    stop("container must be a 'Cel.container'")
  
  n <- length( container )
  
  if (! (refindex %in% 1:n)) stop("invalid reference index for normalization")
  
  refmean <- container[[refindex]]@intensity
  
  for (i in 1:n) {
    container[[i]]@intensity <- normalize.constant(container[[i]]@intensity, refmean)
    container[[i]]@history <- list(name="normalized by constant",
                                 constant=attr(container[[i]]@intensity,"constant"))
    container[[i]]@sd <- NULL
    attr(container[[i]]@intensity, "constant") <- NULL
  }
  return(container)
}       


normalize.constant <- function(data, refmean) {
  thismean <- mean(data)
  r <- data / thismean * refmean
  attr(r,"constant") <- thismean * refmean
}
