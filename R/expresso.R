expresso <- function(afbatch,
                     bg.correct=TRUE,
                     bgcorrect.method = NULL,
                     bgcorrect.param = list(),
                     normalize = TRUE,
                     normalize.method = NULL,
                     normalize.param=list(),
                     pmcorrect.method = NULL,
                     pmcorrect.param = list(),
                     summary.method = NULL,
                     summary.param = list(),
                     summary.subset = NULL,
                     ## ---
                     ##phenodata = NULL,##assume it comes with AffyBatch, samee for MIAME
                     verbose = TRUE,
                     warnings = TRUE,
                     widget = FALSE
                     ) {
  
  
  
  if (widget) {
    require(tkWidgets) || stop("library tkWidgets could not be found !")
  }
  
  nchips <- length(afbatch)

  ###background stuff must be added before normalization!
  
  ## -- background correction method
  if (is.null(bgcorrect.method)) {
    if (widget) {
      bgcorrect.method <- pwidget.selector(bgcorrect.methods,
                                    title = "Method for background correction")
                                    ##,choices.help = NULL)
      #bg.method <- paste("bg.correct", bg.method, sep=".")
    } else {
      stop("bg.method missing")
    }
  }
  
  ## -- normalize.method
  if ((normalize) & (is.null(normalize.method))) {
    if (widget) {
      ##DEBUG: move what's below to 'pwidget.selector'
      n.methods <- normalize.methods(afbatch)
      normalize.method <- pwidget.selector(n.methods,
                                           title = "Method for normalization")
      ##choices.help = paste("normalize", n.methods,  sep="."))##took this out cause netscape popping up was annoying me
      
      rm(n.methods)
    } else {
      stop(paste("normalization method missing. Try one of:",
                 normalize.methods(afbatch), sep="\n"))
    }
  }

  ## -- background correction method
  if (is.null(pmcorrect.method)) {
    if (widget) {
      pmcorrect.method <- pwidget.selector(pmcorrect.methods,
                                    title = "Method for PM correction")
      ##,choices.help = NULL)
      #bg.method <- paste("bg.correct", bg.method, sep=".")
    } else {
      stop("pmcorrect.method missing")
    }
  }
  
  ## -- expression method
  if (is.null(summary.method)) {
    if (widget) {
      helpnames <- paste("generateExprVal.method", generateExprSet.methods, sep="")
      summary.method <- pwidget.selector(generateExprSet.methods,
                                         title = "Method for expression")
                                         ##choices.help = helpnames)
    } else {
      stop("summary.method method missing")
    }
  }
  
  ## -- summary of what will be done
  if (verbose) {
    if (bg.correct){
      cat("background correction:", bgcorrect.method, "\n")
    }
    if (normalize) {
      cat("normalization:", normalize.method, "\n")
    }
    cat("PM/MM correction :", pmcorrect.method, "\n")
    cat("expression values:", summary.method, "\n")
  }

  ##this can be chagnes through sampleNames <-
#   ## --- chip.names
#   if (is.null(chip.names)) {
#     chip.names <- as.character(sapply(as.list(1:nchips), function(x) strsplit(x,"\\.")[[1]][1]))
#   } else {
#     if (length(chip.names) != nchips) {
#       warning("Not the same number of chips than chip names. Assigning names from file.\n")
#       chip.names <- as.character(sapply(as.list(1:nchips), function(x) strsplit(x,"\\.")[[1]][1])) 
#     }
#   }
  
  ## --- reading CDF (if needed)

  ## -- background correct (if needed)
  if (bg.correct) {
    
    if (verbose)
      cat("background correcting...")
    
    afbatch <- do.call("bg.correct", c(alist(afbatch, method=bgcorrect.method), bgcorrect.param))
    
    if (verbose)
      cat("done.\n")
  }
    
  ## -- normalize (if wished)
  if (normalize) {
    
    if (verbose)
      cat("normalizing...")
    
    afbatch <- do.call("normalize",
                       c(alist(afbatch, normalize.method), normalize.param))
    
    if (verbose)
      cat("done.\n")
  }  
  
  eset <- computeExprSet(afbatch, #bg.method=bg.method,
                         summary.method=summary.method, pmcorrect.method= pmcorrect.method,
                         ids=summary.subset,
                         summary.param=summary.param, pmcorrect.param=pmcorrect.param,
                         warnings=warnings)
  
  ##  if (! is.null(phenodata)), ##must assume we get it 
  
  return(eset)
}



