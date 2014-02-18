xy2indices <- function(x, y, nc=NULL, cel=NULL, abatch=NULL, cdf=NULL, xy.offset = NULL) {

  if ( is.null(xy.offset) ) {
    xy.offset <- getOption("BioC")$affy$xy.offset
  }
  
  if (any(x < xy.offset) || any(y < xy.offset))
    stop("Xs and Ys must start at 0 or 1 (please refer to the help file) !")
  
  ct <- sum(c(is.null(nc), is.null(cel), is.null(abatch), is.null(cdf)))
  if (ct != 3)
    stop("One and only one of 'nc', 'cel', 'abatch', 'cdf' should be specified.")
  if (! is.null(cel))
    stop("Cel class no longer supported") #nr <- nrow(intensity(cel))
  if (! is.null((abatch)))
    nc <- ncol(abatch)
  if(!is.null(cdf)){
    require(cdf, character.only = TRUE, quietly = TRUE) ||
            stop(paste(cdf, "package must be installed first"))
    nam <- paste(sub("cdf$", "", cdf), "dim", sep = "")
    nc <- get("NCOL", get(nam))
  }
  
  return( (x - xy.offset) + 1 + nc * (y - xy.offset) )
}


indices2xy <- function(i, nc=NULL, cel=NULL, abatch=NULL, cdf=NULL, xy.offset = NULL) {

  if ( is.null(xy.offset) ) {
    xy.offset <- getOption("BioC")$affy$xy.offset
  }
  
  if (any(i <= 0))
    stop("Indices must start at 0 or 1 (please refer to the help file) !")

  ct <- sum(c(is.null(nc), is.null(cel), is.null(abatch), is.null(cdf)))
  
  if (ct != 3)
    stop("One and only one of 'nc', 'cel', 'abatch', 'cdf' should be specified.")
  if (! is.null(cel))
    stop("Cel class no longer supported")#    nr <- nrow(intensity(cel))
  if (! is.null((abatch)))
    nc <- ncol(abatch)
  if(!is.null(cdf)){
    require(cdf, character.only = TRUE, quietly = TRUE) ||
            stop(paste(cdf, "package must be installed first"))
    nam <- paste(sub("cdf$", "", cdf), "dim", sep = "")
    nc <- get("NCOL", get(nam))
  }
  
  
  ##x <- i %% nr
  ##x[x == 0] <- nr
  ##y <- (i - 1) %/% nr + 1  
  x <- (i  - 1) %% nc + xy.offset
  y <- (i - 1) %/% nc + xy.offset
  xy <- cbind(x, y)
  colnames(xy) <- c("x", "y")
  return(xy)
}
