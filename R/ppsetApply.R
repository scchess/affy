ppsetApply <- function(abatch, FUN, genenames=NULL, ...) {

  if (! inherits(abatch, "AffyBatch"))
    stop("abatch must be inheriting from class AffyBatch")

  if (! inherits(FUN, "function"))
    stop("FUN must be a function")
  
  cdfenv <- getCdfInfo(abatch)

  if (is.null(genenames))
    genenames <- ls(cdfenv)

  ##
  e1 <- new.env(parent = environment(FUN))
  multiassign(names(pData(abatch)), pData(abatch), env = e1)
  environment(FUN) <- e1
  
  ppset <- new("ProbeSet", pm=matrix(), mm=matrix())

  r <- vector("list", length=length(genenames))
  names(r) <- genenames
  
  for (i in seq(along=genenames)) {
    ## use multiget to get NA when genenames[i] not found
    probes.i <- multiget(genenames[i], envir = cdfenv)[[1]]
    if (all(is.na(probes.i)))
      next
    ppset@pm <- intensity(abatch)[probes.i[, 1], , drop=FALSE]
    ppset@mm <- intensity(abatch)[probes.i[, 2], , drop=FALSE]
    ppset@id <- genenames[i]
    r[[i]] <- FUN(ppset, ...)
  }

  return(r)
}

ppsetClusterApply <- function(abatch, FUN, genenames=NULL, ...) {

  if (! inherits(abatch, "AffyBatch"))
    stop("abatch must be inheriting from class AffyBatch")

  if (! inherits(FUN, "function"))
    stop("FUN must be a function")
  
  cdfenv <- getCdfInfo(abatch)

  if (is.null(genenames))
    genenames <- ls(cdfenv)

  ##
  e1 <- new.env(parent = environment(FUN))
  multiassign(names(pData(abatch)), pData(abatch), env = e1)
  environment(FUN) <- e1
  
  ppset <- new("ProbeSet", pm=matrix(), mm=matrix())

  r <- vector("list", length=length(genenames))
  names(r) <- genenames
  
  for (i in seq(along=genenames)) {
    ## use multiget to get NA when genenames[i] not found
    probes.i <- multiget(genenames[i], envir = cdfenv)[[1]]
    if (all(is.na(probes.i)))
      next
    ppset@pm <- intensity(abatch)[probes.i[, 1], , drop=FALSE]
    ppset@mm <- intensity(abatch)[probes.i[, 2], , drop=FALSE]
    ppset@id <- genenames[i]
    r[[i]] <- FUN(ppset, ...)
  }

  return(r)
}


ppset.ttest <- function(ppset, covariate, pmcorrect.fun = pmcorrect.pmonly, ...) {
  probes <- do.call("pmcorrect.fun", list(ppset))
  my.ttest <- function(x) {
    y <- split(x, get(covariate))
    t.test(y[[1]], y[[2]])$p.value
  }
  r <- apply(probes, 1, my.ttest)
  return(r)
}



# make.ppset.logitt <- function(abatch) {
  
#   ppset.logitt <- function(ppset, covariate, pmcorrect.fun = pmcorrect.pmonly, A, N) {
#     probes <- do.call("pmcorrect.fun", list(ppset))
#     probes.logit <- 
#     }

#   return(ppset.logitt)
# }

