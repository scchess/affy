##this function does the expression for our method: RMA. not flexible
rma <- function(object, subset=NULL, verbose=TRUE, normalize=TRUE,
                bg.correct=TRUE, ...){
  
  
  if(bg.correct){
    if(verbose) cat("Background correcting\n")
    ##ADD pm<- to methods!
    object <- bg.correct.rma(object)
  }

  if(normalize){
    if(verbose) cat("Normalizing Data\n")
    object <- normalize(object,"quantiles", pmonly=TRUE)
  }
  
  ## if NULL compute for all
  if (is.null(subset))
    subset <- geneNames(object)

  n <- length(object)
  m <- length(subset)

  c.pps <- new("ProbeSet",pm=matrix(), mm=matrix())
              
  ## matrix to hold expression values
  exp.mat <- matrix(NA, m, n)
  se.mat <- matrix(NA, m, n)
  
  if (verbose) {
    cat(m, "ProbeSets to be processed\n")
    countprogress <- 0
  }
              
  mycall <- as.call(c(getMethod("express.summary.stat", signature=c("ProbeSet","character")),list(c.pps, method="medianpolish",param.method=list())))


  CDFINFO <- getCdfInfo(object) ##do it once!
  for (i in seq(along=subset)) {
                
    id <- subset[i]
    if (verbose) {
      if ( round(m/10) == countprogress) {
        cat(".")
        countprogress <- 0
      }
      else
        countprogress <- countprogress + 1
    }
    loc <- get(id ,envir=CDFINFO)
    l.pm <- loc[, 1]
                
    np <- length(l.pm)
                
    c.pps@pm <- matrix(intensity(object)[l.pm, ],
                       np, n, byrow=TRUE)
                
    mycall[[2]] <- c.pps
    ev <- eval(mycall)
                
    exp.mat[i, ] <- ev$exprs
    se.mat[i,] <- ev$se.exprs
  }

  if (verbose) cat("\n")
  
  dimnames(exp.mat) <- list(subset, sampleNames(object))
  dimnames(se.mat) <- list(subset, sampleNames(object))

  object@exprs <- exp.mat
  object@se.exprs <- se.mat
  class(object) <- "exprSet"
  return(object)
}
  





