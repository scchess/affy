normalize.Cel.container.quantiles <- function(listcel) {
  ## make the matrix (not very memory friendly... will be better when HDF5,
  ## or an another mecanism is around...)
  ## note to Rafael: I hope I did not butcher your code too much... this is just
  ## an example on how to get the function to the new objects scheme
  ## (note on the note: you may want to create a 'Plob' object the very same way I
  ## did with 'Cel' and the corresponding generic method (see 'Cel.R'), and of course
  ## the 'normalize.Plob.yyy' methods (see 'normalize.invariantset.R')).
  ##
  ## Laurent
  
  cols <- length(listcel)
  rows <- length(intensity(listcel[[1]])) # assuming all the Cel are of the same size in listcel
  chipdim <- dim(intensity(listcel[[1]]))
  
  ## this is may be too much.. I am taking every intensity in there... something like
  ## taking only the PM (see invariantset) could be what you want (this is only a trial,
  ## know better your own algorithm that I will ever... modify if needed).  
  x <- matrix(0,rows,cols)
  for (i in 1:cols) x[,i] <- c(intensity(listcel[[i]]))

  x <- normalize.quantiles(x)

  for (i in 1:cols) {
    intensity(listcel[[i]])  <- matrix(x[,i],chipdim[1], chipdim[2])
    history(listcel[[i]])$name <- "normalized by quantiles"
    spotsd(listcel[[i]]) <- matrix() # set 'sd' to nothing (meaningless after normalization)
  }

  return(listcel)
}

normalize.Plob.quantiles <- function(object,...) {

  x <- normalize.quantiles(rbind(pm(object),mm(object)))
  n <- dim(x)[1]/2
  pm(object) <- x[1:n,] 
  mm(object) <- x[(n+1):(2*n),]
  return(object)

}

normalize.quantiles <- function(x){

  rows <- dim(x)[1]
  cols <- dim(x)[2]
  
  matrix(.C("qnorm_c",as.double(as.vector(x)),as.integer(rows),as.integer(cols))[[1]],rows,cols)
}
