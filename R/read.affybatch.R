read.affybatch <- function(..., filenames=character(0),
                           ##sd=FALSE,
                           phenoData=new("phenoData"),
                           description=NULL,
                           notes="",
                           compress=getOption("BioC")$affy$compress.cel,
                           rm.mask=FALSE, rm.outliers=FALSE, rm.extra=FALSE,
                           hdf5=FALSE, hdf5FilePath=NULL,
                           widget = FALSE,
                           verbose=FALSE) {
  
  auxnames <- as.list(substitute(list(...)))[-1]
  widgetfiles <- character(0)
  if (widget) {
    widgetfiles <- fileBrowser(textToShow="Choose CEL files",
                               testFun=hasSuffix("[cC][eE][lL]"))
  }
  filenames <- .Primitive("c")(filenames, auxnames, widgetfiles)

  if (widget && is.null(description)) {
    description <- tkMIAME()
  }
  if (is.null(description))
    description <- new("MIAME")
  
  n <- length(filenames)
  
  ## error if no file name !
  if (n == 0)
    stop("No file name given !")
  
  ## read the first file to see what we have
  if (verbose) cat(1, "reading",filenames[[1]],"...")
  
  cel <- read.celfile(filenames[[1]],
                      ##sd=sd,
                      compress=compress,
                      rm.mask = rm.mask,
                      rm.outliers = rm.outliers,
                      rm.extra = rm.extra)
  if (verbose) cat("done.\n")
  
  dim.intensity <- dim(intensity(cel))
  ##if (sd)
  ##  dim.sd <- dim.intensity
  ##else
  ##  dim.sd <- c(1,1)
  
  if (hdf5) {
    require(rhdf5) || stop("The package rhdf5 is required !")
    if (is.null(hdf5FilePath))
      stop("A path for tmp files must be specified")
    if (! is.na(file.info(hdf5FilePath)$size)) {
      warning(paste("The file \"", hdf5FilePath, "\" already exists !", sep=""))
    }
    conty <- new.AffyBatch.hdf5(n, dim.intensity,
                                hdfile.group="raw",
                                hdfile.name=hdf5FilePath,
                                cdfName = cel@cdfName
                                )
    
  } else {
    conty <- new("AffyBatch",
                 intensity  = array(NA, dim=c(dim.intensity, n)),
                 ##sd = array(NA, dim=dim.sd),
                 cdfName    = cel@cdfName,
                 phenoData  = phenoData,
                 nrow       = dim.intensity[1],
                 ncol       = dim.intensity[2],
                 nexp       = n,
                 description= description,
                 notes      = notes,
                 history    = vector("list", length=n))
    dimnames(intensity(conty)) <- list(NULL, NULL, rep("", n))
  }
  
  intensity(conty)[, , 1] <- intensity(cel)
  
  ##if (sd)
  ##  spotsd(conty)[, , 1] <- spotsd(cel)
  
  c.names <- rep("", n)
  c.names[1] <- cel@name
  ##outliers(conty)[[1]] <- outliers(cel)
  ##masks(conty)[[1]] <- masks(cel)
  history(conty)[[1]] <- history(cel)

  ## finish if only one file
  if (n == 1)
    return(conty)

  for (i in 2:n) {
    
    if (verbose) cat(i, "reading",filenames[[i]],"...")
    cel <- read.celfile(filenames[[i]],
                        ##sd=sd,
                        compress=compress, rm.mask=rm.mask,
                        rm.outliers=rm.outliers, rm.extra=rm.extra)
    if (dim(intensity(cel)) != dim.intensity)
      stop(paste("CEL file dimension mismatch !\n(file",filenames[[i]],")"))
    if (verbose) cat("done.\n")

    if (cel@cdfName != conty@cdfName)
      warning(paste("cdfName mismatch !\n(", filenames[[i]], ")"))
    
    intensity(conty)[, , i] <- intensity(cel)
    dimnames(intensity(conty))[[3]][i] <- cel@name
    ##if (sd)
    ##  spotsd(conty)[, , i] <- spotsd(cel)
    c.names[i] <- cel@name
    ##outliers(conty)[[i]] <- outliers(cel)
    ##masks(conty)[[i]] <- masks(cel)
    history(conty)[[i]] <- history(cel)
  }
  dim(intensity(conty)) <- c(prod(dim.intensity), n)
  chipNames(conty) <- c.names 

  return(conty)
  
}
