expresso <- function(CDFfile = NULL,
                     CELfiles = NULL,
                     compress.cdf = FALSE,
                     compress.cel = FALSE,
                     rm.mask = FALSE,
                     rm.outliers = FALSE,
                     rm.extra = FALSE,
                     bg.method=NULL,
                     chip.names = NULL,
                     phenodata=NULL,
                     annotation=NULL,
                     notes=NULL,
                     description=NULL,
                     subset=NULL,
                     summary.stat=medianpolish,
                     normalize=T,
                     normalize.method="quantiles",
                     eset.method=NULL,
                     verbose = T,
                     widget=F,
                     ...) {

  require(rhdf5) || stop("library rhdf5 could not be found !")
  
  hdf5FilePath <- "test.hd5"
  
  if (widget) {
    require(tkWidgets) || stop("library tkWidgets could not be found !")
  }

  ## --- CDF
  if (is.null(CDFfile)) {
    if (widget) {
      CDFfile <- fileBrowser(textToShow="Choose one CDF file",nSelect=1,
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
      normalize.method <- pwidget.normalize(new("Cel.container"))
    } else {
      stop("normalization method missing")
    }
  }

  ## -- summary of what will be done
  if (verbose) {
    if (normalize) {
      cat("normalization:", normalize.method, "\n")
    }
    cat("background correction:", bg.method, "\n")
    cat("expression values:", eset.method, "\n")
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

  ## --- reading CDF
  if (verbose) cat("reading CDF file...")
  cdf <- read.cdffile(CDFfile, compress = compress.cdf)
  if (verbose) cat("done.\n")

  ## --- reading CELs
  if (verbose) cat("reading ", length(CELfiles), "CEL file(s):\n")
  listcel <- read.container.celfile(filenames=CELfiles,
                                    compress=compress.cel,
                                    rm.mask=rm.mask,
                                    rm.outliers=rm.outliers,
                                    rm.extra=rm.extra,
                                    sd=FALSE,
                                    hdf5=TRUE,
                                    hdf5FilePath=hdf5FilePath,
                                    verbose=verbose)

  ## -- normalize (if wished)
  if (normalize) {
    
    if (verbose) cat("normalizing...")
    ## trick: make a copy of listcel (as it will be modifed by the
    ## normalization routines)
    listcel.n <- new.Cel.container.hdf5(n,
                                        hdf5FilePath,
                                        dim.intensity=dim(intensity(listcel[[1]])),
                                        hdfile.group="normalized")
    listcel.n@name <- listcel@name
    outliers(listcel.n) <- outliers(listcel)
    masks(listcel.n) <- masks(listcel)
    history(listcel.n) <- history(listcel)
    ## DEBUG: copy the contents 
    intensity(listcel.n)[, , ] <- intensity(listcel)[, , ]
    if (sd)
      spotsd(listcel.n)[, , ] <- spotsd(listcel)[, , ]
    else
      spotsd(listcel.n)[, ] <- spotsd(listcel)[, ]
    
    listcel.n <- normalize(listcel.n, cdf, method=normalize.methods, hdf5=TRUE)
    
    listcel <- listcel.n
    
    if (verbose) cat("done.\n")
  }
  
  eset <- generateExprSet(listcel, cdf, method=eset.method, bg.correct=bg.method)

  return(eset)
}
