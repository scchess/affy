expresso <- function(cdf = NULL,
                     CELfiles = NULL,
                     compress.cdf = FALSE,
                     compress.cel = FALSE,
                     rm.mask = FALSE,
                     rm.outliers = FALSE,
                     rm.extra = FALSE,
                     ## ---
                     normalize = TRUE,
                     normalize.method = NULL,
                     normalize.param=list(),
                     bg.method = NULL,
                     bg.param = list(),
                     summary.method = NULL,
                     summary.param = list(),
                     summary.subset = NULL,
                     chip.names = NULL,
                     ## ---
                     phenodata = NULL,
                     verbose = T,
                     widget = F,
                     hdf5 = F) {
  

  ## --- temp file (if hdf5)
  if (hdf5) {
    require(rhdf5) || stop("library rhdf5 could not be found !")
    
    hdf5FilePath <- .getTmpFileName(".hf5")
    
    if (verbose)
      cat("temporary hdf5 file:",hdf5FilePath,"\n")
  } else {
    hdf5FilePath <- "dummy"
  }
  
  if (widget) {
    require(tkWidgets) || stop("library tkWidgets could not be found !")
  }

  ## --- CDF
  if (is.null(cdf)) {
    if (widget) {
      cdf <- fileBrowser(textToShow="Choose one CDF file",nSelect=1,
                         testFun=hasSuffix("[cC][dD][fF]"))
    } else {
      stop("CDF file missing")
    }
  }

  ## --- CEL
  if (is.null(CELfiles)) {
    if (widget) {
      CELfiles <- fileBrowser(textToShow="Choose CEL files",
                              testFun=hasSuffix("[cC][eE][lL]"))
    } else {
      stop("CEL file missing")
    }
  }
  nchips <- length(CELfiles)

  ## -- normalize.method
  if ((normalize) & (is.null(normalize.method))) {
    if (widget) {
      ##DEBUG: move what's below to 'pwidget.selector'
      n.methods <- normalize.methods(new("Cel.container"))
      normalize.method <- pwidget.selector(n.methods,
                                           title = "Method for normalization",
                                           choices.help = paste("normalize", n.methods, sep="."))
      rm(n.methods)
    } else {
      stop("normalization method missing")
    }
  }

  ## -- background correction method
  if (is.null(bg.method)) {
    if (widget) {
      bg.method <- pwidget.selector(bg.correct.methods,
                                    title = "Method for background correction",
                                    choices.help = NULL)
      #bg.method <- paste("bg.correct", bg.method, sep=".")
    } else {
      stop("bg.method missing")
    }
  }
  
  ## -- expression method
  if (is.null(summary.method)) {
    if (widget) {
      helpnames <- paste("generateExprVal.method", generateExprSet.methods, sep="")
      summary.method <- pwidget.selector(generateExprSet.methods,
                                         title = "Method for expression",
                                         choices.help = helpnames)
    } else {
      stop("summary.method method missing")
    }
    
  }
  ## -- summary of what will be done
  if (verbose) {
    if (normalize) {
      cat("normalization:", normalize.method, "\n")
    }
    cat("background correction:", bg.method, "\n")
    cat("expression values:", summary.method, "\n")
  }
  
  ## --- chip.names
  if (is.null(chip.names)) {
    chip.names <- as.character(sapply(CELfiles, function(x) strsplit(x,"\\.")[[1]][1]))
  } else {
    if (length(chip.names) != nchips) {
      warning("Not the same number of chips than chip names. Assigning names from file.\n")
      chip.names <- as.character(sapply(CELfiles, function(x) strsplit(x,"\\.")[[1]][1])) 
    }
  }

  ## --- reading CDF (if needed)
  if (class(cdf) != "Cdf") {
    if (verbose) cat("reading CDF file...")
    
    ##cdf <- try(read.cdffile(CDFfile, compress = compress.cdf))
    ##if (inherits(cdf,"try-error")) {
    ##  if (verbose) cat("(trying again with/without compression)...")
    ##  cdf <- try(read.cdffile(CDFfile, compress = !compress.cdf))
    ##}
    cdf <- read.cdffile(cdf, compress = compress.cdf)
    
    if (verbose) cat("done.\n")
  }
  
  ## --- reading CELs
  if (verbose) cat("reading ", length(CELfiles), "CEL file(s):\n")

  if (hdf5) {
    on.exit(cat("unlinking temporary file", hdf5FilePath, ".\n"), add=TRUE)
    on.exit(unlink(hdf5FilePath), add=TRUE)
  }
  
#   listcel <- try(read.container.celfile(filenames=CELfiles,
#                                         compress=compress.cel,
#                                         rm.mask=rm.mask,
#                                         rm.outliers=rm.outliers,
#                                         rm.extra=rm.extra,
#                                         sd=FALSE,
#                                         hdf5=TRUE,
#                                         hdf5FilePath=hdf5FilePath,
#                                         verbose=verbose)
#                  )
#   if (inherits(listcel,"try-error")) {
#     if (verbose) cat("(trying again with/without compression)...")
#     listcel <- try(read.container.celfile(filenames=CELfiles,
#                                           compress=compress.cel,
#                                           rm.mask=rm.mask,
#                                           rm.outliers=rm.outliers,
#                                           rm.extra=rm.extra,
#                                           sd=FALSE,
#                                           hdf5=TRUE,
#                                           hdf5FilePath=hdf5FilePath,
#                                           verbose=verbose)
#                    )
#   }
  listcel <- read.container.celfile(filenames=CELfiles,
                                    compress=compress.cel,
                                    rm.mask=rm.mask,
                                    rm.outliers=rm.outliers,
                                    rm.extra=rm.extra,
                                    sd=FALSE,
                                    hdf5=hdf5,
                                    hdf5FilePath=hdf5FilePath,
                                    verbose=verbose)
  
  
  ## -- normalize (if wished)
  if (normalize) {

    if (hdf5) {
      hdf5FilePath2 <- .getTmpFileName(".hd5")
      on.exit(cat("unlinking temporary file", hdf5FilePath2, ".\n"), add=TRUE)
      on.exit(unlink(hdf5FilePath2), add=TRUE)
      
      if (verbose) {
        cat("temporary hdf5 file:",hdf5FilePath2,"\n")
        cat("making a copy of the raw values...")
      }
      
      ## trick: make a copy of listcel (as it will be modifed by the
      ## normalization routines)
      
      listcel.n <- convert2hdf5.Cel.container(listcel,
                                              hdf5FilePath2,
                                              hdfile.group = "normalized")
      
      if (verbose) {
        cat("done.\n")
        cat("normalizing...")
      }
      
      listcel.n <- do.call("normalize", c(list(listcel.n, cdf, method=normalize.method), normalize.param))
      
      listcel <- listcel.n
      
    } else {
      ## hdf5 == FALSE
      if (verbose)
        cat("normalizing...")
      listcel <- do.call("normalize", c(list(listcel, cdf, method=normalize.method), normalize.param))
      if (verbose)
        cat("done.\n")
    }
    
    if (verbose) cat("done.\n")
  }
  
  eset <- generateExprSet(listcel, cdf, method=summary.method, bg.correct=bg.method, ids=summary.subset)

  if (! is.null(phenodata))
    phenoData(eset) <- phenodata
  
  return(eset)
}
