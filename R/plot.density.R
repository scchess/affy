plotDensity <- function(mat,
                        ylab="density", xlab="x", ...) {
  
  x.density <- apply(mat, 2, density)

  all.x <- do.call("cbind", lapply(x.density, function(x) x$x))
  all.y <- do.call("cbind", lapply(x.density, function(x) x$y))
  
  matplot(all.x, all.y, ylab=ylab, xlab=xlab, ...)
}
 

plotDensity.AffyBatch <- function(x, col=rainbow(length(x)), log=TRUE,
                                  which=c("pm","mm","both"),
                                  ylab="density",
                                  xlab=NULL,
                                  ...){
  
  Index <- unlist(indexProbes(x, which=which))
  
  x <- intensity(x)[Index, ]
  
  if(log){
    x <- log2(x)
    if(is.null(xlab)) xlab <- "log intensity"
  }
  else  if(is.null(xlab)) xlab <- "intensity"
  
  
  plotDensity(x, ylab=ylab, xlab=xlab, ...)
}
