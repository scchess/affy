## (see normalize.invariantset.R to know more....)

normalize.Cel.constant <- function(listcel, refindex=1) {
  if (! inherits(listcel, "Cel.container"))
    stop("listcel must be a 'Cel.container'")
  
  n <- length( listcel )
  
  if (! (refindex %in% 1:n)) stop("invalid reference index for normalization")
  
  refmean <- listcel[[refindex]]@intensity
  
  for (i in 1:n) {
    listcel[[i]]@intensity <- normalize.constant(listcel[[i]]@intensity, refmean)
    listcel[[i]]@history <- list(name="normalized by constant",
                                 constant=attr(listcel[[i]]@intensity,"constant"))
    listcel[[i]]@sd <- NULL
    attr(listcel[[i]]@intensity, "constant") <- NULL
  }
  return(listcel)
}       


normalize.constant <- function(data, refmean) {
  thismean <- mean(data)
  r <- data / thismean * refmean
  attr(r,"constant") <- thismean * refmean
}
