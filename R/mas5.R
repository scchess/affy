mas5 <- function(object,normalize=TRUE,sc = 500, analysis = "absolute",...){
  res <- expresso(object,bgcorrect.method="mas",pmcorrect.method="mas",normalize=FALSE,summary.method="mas",...) 
  if(normalize) 
    res <- affy.scalevalue.exprSet(res,sc=sc,analysis=analysis)
  return(res)
}

mas5calls.ProbeSet <- function(object, 
                               tau=0.015, alpha1=0.04, alpha2=0.06,
                               exact.pvals=FALSE, cont.correct=FALSE) {
  pms <- pm(object)
  mms <- mm(object)
  calls <- vector("character",ncol(pms))
  pvals <- vector("numeric",ncol(pms))
  for(i in 1:ncol(pms)){
    mat <- cbind(pms[,i],mms[,i])
    res <- mas5.detection(mat,tau=tau, alpha1=alpha1, alpha2=alpha2,
                           exact.pvals=exact.pvals, cont.correct=cont.correct)
    calls[i] <- res$call
    pvals[i] <- res$pval
  }
  return(list(call=calls,pval=pvals))
}

mas5calls.AffyBatch <- function(object, ids=NULL, verbose=TRUE,
                                tau=0.015, alpha1=0.04, alpha2=0.06,
                                exact.pvals=FALSE, cont.correct=FALSE) {
  
  n <- length(object)
  
  ## if NULL compute for all
  if (is.null(ids))
    ids <- geneNames(object)
  
  m <- length(ids)
  
  exp.mat <- matrix(NA, m, n)
  se.mat <- matrix(NA, m, n)

  if (verbose) {
    cat(m, "ids to be processed\n")
    countprogress <- 0
  }

  ## loop over the ids
  
  CDFINFO <- getCdfInfo(object) ##do it once!

  for (i in seq(along=ids)) {

    id <- ids[i]

    if (verbose) {
      if ( round(m/10) == countprogress) {
        cat(".")
        countprogress <- 0
      }
      else
        countprogress <- countprogress + 1
    }

    ## locations for an id
    loc <- get(id, envir=CDFINFO)
    l.pm <- loc[, 1]
    if (ncol(loc) == 2)
      l.mm <- loc[ ,2]
    else
      stop("No MMs available for id:",id,"\n")

    pms <- intensity(object)[l.pm, ,drop=FALSE]
    mms <- intensity(object)[l.mm, ,drop=FALSE]
    for(j in 1:n){
      mat <- cbind(pms[,j],mms[,j])

      res <- mas5.detection(mat, tau=tau, alpha1=alpha1, alpha2=alpha2,
                           exact.pvals=exact.pvals, cont.correct=cont.correct)
      exp.mat[i,j] <- res$call
      se.mat[i,j] <- res$pval
    }
  }
  if (verbose) cat("\n")
  
  ## instance exprSet
  dimnames(exp.mat) <- list(ids, sampleNames(object))
  dimnames(se.mat) <- list(ids, sampleNames(object))
  eset <- new("exprSet",
              exprs=exp.mat,
              se.exprs=se.mat,
              phenoData=phenoData(object),
              description=description(object),
              annotation=annotation(object),
              notes=c(notes(object)))
  return(eset)
}


mas5.detection <- function(mat, tau=0.015, alpha1=0.04, alpha2=0.06,
                           exact.pvals=FALSE, cont.correct=FALSE) { 

  mat.r <- (mat[,1]-mat[,2])/(mat[,1]+mat[,2])
  ## CONSTANTS
  saturation.point <- 46000			# not a user parameter
     
     ## SANITY CHECKING
    if ( !is.matrix(mat) || length(dim(mat))!=2 || dim(mat)[2]!=2 ||
         dim(mat)[1] < 1 || !is.numeric(mat) )
        stop("Invalid mat matrix.")
    if ( !is.numeric(tau) )
        stop("Invalid tau.")
    if ( !is.numeric(alpha1) || !is.numeric(alpha2) ||
          alpha1 <= 0 || alpha1 >= alpha2 || alpha2 >= 0.5 )
        stop("Invalid alpha1 or alpha2.")
    if ( !is.logical(exact.pvals) )
        stop("Invalid exact.pvals.")
    if ( !is.logical(cont.correct) )
        stop("Invalid cont.correct.")
        
     ## DEALING WITH SATURATION; COMPUTING THE P-VALUE
     ## According to the Bioinformatics paper:
     ## * If all MM's are saturated, then call present
     ## * Otherwise discard pairs with a saturated MM
     ## According to the Affymetrix whitepaper:
     ## * If all probe-pairs are saturated, then call present with pval=0
     ## * If an MM is saturated, then we discard the pair
     ## * If a PM and MM are within tau of each other, we discard the pair
     ## So we're going with:
     ## * If all MM's are saturated, set pval=0 and don't use Wilcoxon
     ## * Discard probe-pairs when MM is saturated or the PM,MM are within tau
     ##   of each other
     ## * Compute the p-value using Wilcoxon's signed rank test on the retained
     ##   probe-pairs
    is.mm.saturated <- function(probe.pair, saturation.point)
        probe.pair[2] >= saturation.point
    is.retained <- function(probe.pair, saturation.point, tau)
        !(is.mm.saturated(probe.pair,saturation.point) ||
          abs(diff(probe.pair)) <= tau)
    if ( all(apply(mat,1,is.mm.saturated,saturation.point)) )
        pval <- 0
    else {
      retained <- apply(mat, 1, is.retained, saturation.point, tau)
      pval <- wilcox.test(mat.r[retained],
                          alternative="greater", mu=tau, paired=FALSE,
                          exact=exact.pvals, correct=cont.correct,
                          conf.int=FALSE)$p.value
    }
  
     ## DETECTION CALL
    if ( pval < 0 || pval > 1 )
        warning("Computed an unusual p-value outside the range [0,1].")
    if ( pval < alpha1 )
        call <- "P"
    else if ( pval < alpha2 )
        call <- "M"
    else
        call <- "A"
    
     ## DONE
    return(list(pval=pval, call=call))
}

