expresso <- function(afbatch,
                     normalize = TRUE,
                     normalize.method = NULL,
                     normalize.param=list(),
                     bg.correct=TRUE,
                     bg.method = NULL,
                     ##bg.param = list(),
                     summary.method = NULL,
                     summary.param = list(),
                     summary.subset = NULL,
                     ## ---
                     ##phenodata = NULL,##assume it comes with AffyBatch, samee for MIAME
                     verbose = TRUE,
                     widget = FALSE
                     ) {
  
  
  
  if (widget) {
    require(tkWidgets) || stop("library tkWidgets could not be found !")
  }
  
  nchips <- length(afbatch)

  ###background stuff must be added before normalization!
  
  ## -- background correction method
  if (is.null(bg.method)) {
    if (widget) {
      bg.method <- pwidget.selector(bg.correct.methods,
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
                                           ###choices.help = paste("normalize", n.methods,  sep="."))##took this out cause netscape popping up was annoying me
                                          
    rm(n.methods)
    } else {
      stop(paste("normalization method missing. Try one of:",
                 normalize.methods(afbatch), sep="\n"))
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
      cat("background correction:", bg.method, "\n")
    }
    if (normalize) {
      cat("normalization:", normalize.method, "\n")
    }
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

  ## -- background correcet (if needed)
  if (bg.correct) {
    
    if (verbose)
      cat("background correcting...")
    afbatch <- do.call(bg.method, c(list(afbatch)))
    if (verbose)
      cat("done.\n")
  }
    
  ## -- normalize (if wished)
  if (normalize) {
    
    if (verbose)
      cat("normalizing...")
    afbatch <- do.call("normalize", c(list(afbatch, normalize.method), normalize.param))
    if (verbose)
      cat("done.\n")
  }
  
  
  eset <- computeExprSet(afbatch, #bg.method=bg.method,
                         summary.method=summary.method, ids=summary.subset)
  
  ##  if (! is.null(phenodata)), ##must assume we get it 

  
  return(eset)
}



