xy2indices <- function(x, y, nc=NULL, cel=NULL, abatch=NULL) {
  if (any(x <= 0) || any(y <= 0))
    stop("Xs and Ys must start at 1 (please refer to the help file) !")
  ct <- sum(c(is.null(nc), is.null(cel), is.null(abatch)))
  if (ct != 2)
    stop("One and only one of 'nc', 'cel', 'abatch' should be specified.")
  if (! is.null(cel))
    nc <- ncol(intensity(cel))
  if (! is.null((abatch)))
    nc <- ncol(abatch)
  
  return(x + nc * (y - 1))
}

indices2xy <- function(i, nr=NULL, cel=NULL, abatch=NULL) {
  if (any(i)<= 0)
    stop("Indices must start at 1 (please refer to the help file) !")
  
  ct <- sum(c(is.null(nr), is.null(cel), is.null(abatch)))
  
  if (ct != 2)
    stop("One and only one of 'nc', 'cel', 'abatch' should be specified.")
  if (! is.null(cel))
    nr <- nrow(intensity(cel))
  if (! is.null((abatch)))
    nr <- nrow(abatch)
  
  x <- i %% nr
  x[x == 0] <- nr
  y <- i %/% nr + 1
  xy <- cbind(x, y)
  colnames(xy) <- c("x", "y")
  return(xy)
}
