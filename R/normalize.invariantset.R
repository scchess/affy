## Laurent 2002

## The idea would be to have the function related to one particular normalization
## technique into one single file (to avoid to run after many different files).
##
## Strongly suggested naming convention:
## (This is not that particular naming convention that matters but that there is A naming
## convention... reasons upon requests for the skeptics...).
## All the functions in such a 'scaling function' file must be like
## 'normalize.xxx.yyy' or 'normalize.yyy'. 'xxx' has to be whether 'Cel' for Affymetrix worshipers,
## whether 'cdna' cDNA arrays devotees (other names will come if other popular formats show up).
## NOTE TO SANDRINE: You are the cDNA array objects guru... it's up to you there... comments ?
##'yyy' has to be the name of the scaling function.

## As you probably already suspect it from above, there is one single function 'normalize.yyy'
## per file. This last function is the data-format-independant algorithm.
## The 'normalize.xxx.yyy' function will communicate with the
## 'normalize.yyy' function according the nature of their data. Of course it may exist functions
## that are not suitable for a particular format 'xxx'. The corresponding function will just
## be absent.
## note: In the case of only one format 'xxx' is by a particular 'yyy' method,
## it may still be wise to have a 'normalize.yyy' function (1- for the future possible formats to
## come, 2- for the re-usability of the code).
##
## The attribute 'history' of the 'Cel' object is used to store informations about the
## procedure (note: they come from the 'scale.yyy' as an attribute). This is probably
## a temporary (but decent) solution.
##
## Careful readers will notice that while the function 'normalize.xxx.yyy' works on a container
## of all the cel objects to scale, the 'normalize.yyy' works on a pair of dataset, one being the
## set to scale and the other the being the reference. This is an attempt to minimize memory
## consumption as this algorithm works sequencially on the data to scale (more careful readers
## may object that this may not change much in the memory usage, to whom I would reply that well
## may be but it would let one make a routine to read sequencially the files and normalize them
## more easily this... anyways HDF5 is coming so those are peculiar considerations... but in the
## meanwhile I can scale 50 chips on my computer at home...=) ).
##
## See Normalize.Cel.Rd to know more about how will one scale painlessly...




## DEBUG
## requires package modreg
library(modreg)
## needed ?


normalize.Plob.invariantset <- function(container, prd.td=c(0.003,0.007), progress=FALSE) {

  nc <- ncol(pm(container))
  np <- nrow(pm(container))
  
  ## take as a reference the array having the median overall intensity
  l <- ncol(pm(container))
  m <- vector("numeric", length=l)
  for (i in 1:l)
    m[i] <- mean(pm(container)[,i])
  refindex <- trunc(median(rank(m)))
  rm(m,l)
  
  ## loop over the CEL files and normalize them
  for (i in (1:nc)[-refindex]) {

    if (progress) cat("normalizing array", attr(container[[i]], "name"), "...")
    
    ##temporary
    tmp <- normalize.invariantset(pm(container)[, i],
                                  pm(container)[, refindex],
                                  prd.td)
                         
    
    # pm first...
    pm(container)[, i] <- as.numeric(approx(tmp$n.curve$y, tmp$n.curve$x,
                                          xout=container@pm[,i])$y, rule=2)
    # then mm... (note: I am not quite sure whether MMs should be 'normalized' or discarded...)
    mm(container)[, i] <- as.numeric(approx(tmp$n.curve$y, tmp$n.curve$x,
                                           xout=mm(container)[,i])$y, rule=2)
    
    ##container@notes <- list(name="normalized by invariant set", invariantset=tmp$i.set)
    container@notes <- "normalized by invariant set"
    
    if (progress) cat("done.\n")
    
  }
  return(container)
}

normalize.Cel.container.invariantset <- function(cel.container, f.cdf, prd.td=c(0.003,0.007), progress=FALSE) {
  
  if (! inherits(cel.container, c("Cel.container","Cel.container.hdf5")))
    stop("cel.container must inherits from 'Cel.container.abstract'")

  if(is.null(f.cdf))
    stop("You need to specify a Cdf object")
  ##DEBUG
  ## test refindex

  i.pm <- pmormm(f.cdf)                               # boolean to find the PM probes
  i.pm[is.na(i.pm)] <- FALSE                          # mark as FALSE the unknown probes too

  np <- sum(i.pm)                                     # number of PM probes
  nc  <-  length(cel.container)                       # number of CEL files
  
  # take as a reference the array having the median overall intensity
  m <- vector("numeric", length=nc)
  for (i in 1:nc)
    m[i] <- mean(intensity(cel.container)[, , i][i.pm])
  refindex <- trunc(median(rank(m)))
  rm(m)           

  if (progress) cat("Data from",cel.container@name[refindex],"used as baseline.\n")

  set.na.spotsd(cel.container)
  
  ## loop over the CEL files and normalize them
  for (i in (1:nc)[-refindex]) {
  
    if (progress) cat("normalizing array", cel.container@name[i], "...")

    mydim <- dim(intensity(cel.container)[, , i])
    
    ##temporary
    tmp <- normalize.invariantset(c(intensity(cel.container)[, , i])[i.pm],
                                  c(intensity(cel.container)[, , refindex])[i.pm],
                                  prd.td)
    i.set <- which(i.pm)[tmp$i.set]
    tmp <- array(as.numeric(approx(tmp$n.curve$y, tmp$n.curve$x,
                                   xout=intensity(cel.container)[, , i], rule=2)$y),
                 mydim)
    attr(tmp,"invariant.set") <- NULL
    intensity(cel.container)[, , i] <- tmp

    ## storing information about what has been done
    history(cel.container)[[i]] <- list(name="normalized by invariant set",
                                       invariantset=i.set)
    
    if (progress) cat("done.\n")
    
  }
  history(cel.container)[[refindex]] <- list(name="reference for the invariant set")
  
  return(cel.container)
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







