read.affybatch <- function(..., filenames=character(0),
                           ##sd=FALSE,
                           phenoData=new("phenoData"),
                           description=NULL,
                           notes="",
                           compress = getOption("BioC")$affy$compress.cel,
                           rm.mask = FALSE, rm.outliers=FALSE, rm.extra=FALSE,
                           verbose = FALSE) {
  
  auxnames <- as.list(substitute(list(...)))[-1]
  filenames <- .Primitive("c")(filenames, auxnames)
  
  n <- length(filenames)
  
  ## error if no file name !
  if (n == 0)
    stop("No file name given !")
  
  pdata <- pData(phenoData)
  ##try to read sample names form phenoData. if not there use CEL filenames
  if(dim(pdata)[1] != n) {
    ##if empty pdata filename are samplenames
    warning("Incompatible phenoData object. Created a new one.\n")
    
    samplenames <- sub("^/?([^/]*/)*", "", unlist(filenames), extended=TRUE)
    pdata <- data.frame(sample=1:n, row.names=samplenames)
    phenoData <- new("phenoData",pData=pdata,varLabels=list(sample="arbitrary numbering"))
  }
  else samplenames <- rownames(pdata)
  
  if (is.null(description))
    {
      description <- new("MIAME")
      description@preprocessing$filenames <- filenames
      description@preprocessing$affyversion <- library(help=affy)$info[[2]][[2]][2]
    }
  ## read the first file to see what we have
  if (verbose) cat(1, "reading",filenames[[1]],"...")
  
  cel <- read.celfile(filenames[[1]],
                      ##sd=sd,
                      compress=compress,
                      rm.mask = rm.mask,
                      rm.outliers = rm.outliers,
                      rm.extra = rm.extra)
  if (verbose) cat("done.\n")
  
  ##now we use the length
  dim.intensity <- dim(intensity(cel))
  ##and the cdfname as ref
  ref.cdfName <- cel@cdfName
  
  ##if (sd)
  ##  dim.sd <- dim.intensity
  ##else
  ##  dim.sd <- c(1,1)

  if (verbose)
    cat(paste("instanciating an AffyBatch (intensity a ", prod(dim.intensity), "x", length(filenames), " matrix)...", sep=""))


  conty <- new("AffyBatch",
               exprs  = array(NaN, dim=c(prod(dim.intensity), n), dimnames=list(NULL, samplenames)),
               ##se.exprs = array(NaN, dim=dim.sd),
               cdfName    = cel@cdfName,
               phenoData  = phenoData,
               nrow       = dim.intensity[1],
               ncol       = dim.intensity[2],
               annotation = cleancdfname(ref.cdfName, addcdf=FALSE),
               description= description,
               notes      = notes)
  ##           history    = vector("list", length=n)) we need to put this in MIAME
  
  if (verbose)
    cat("done.\n")

  # intensity(conty)[, 1] <- c(intensity(cel))
  ##if (sd)
  ##  spotsd(conty)[, , 1] <- spotsd(cel)
  
  ##outliers(conty)[[1]] <- outliers(cel)
  ##masks(conty)[[1]] <- masks(cel)
  ##history(conty)[[1]] <- history(cel) ###this must be done through MIAME
  ival <- intensity(conty)
  ival[, 1] <- c(intensity(cel))
  
  for (i in (1:n)[-1]) {
    
    if (verbose) cat(i, "reading",filenames[[i]],"...")
    cel <- read.celfile(filenames[[i]],
                        ##sd=sd,
                        compress=compress, rm.mask=rm.mask,
                        rm.outliers=rm.outliers, rm.extra=rm.extra)
    
    if (any(dim(intensity(cel)) != dim.intensity))
      stop(paste("CEL file dimension mismatch !\n(file",filenames[[i]],")"))
    if (verbose) cat("done.\n")
    
    if (cel@cdfName != ref.cdfName)
      warning(paste("\n***\nDetected a mismatch of the cdfName: found ", cel@cdfName,
                    ", expected ", ref.cdfName, "\nin file number ", i, " (", filenames[[i]], ")\n",
                    "Please make sure all cel files belong to the same chip type!\n***\n", sep=""))
    
    #intensity(conty)[, i] <- c(intensity(cel))
    
    ##if (sd)
    ##  spotsd(conty)[, , i] <- spotsd(cel)
    
    ##outliers(conty)[[i]] <- outliers(cel)
    ##masks(conty)[[i]] <- masks(cel)
    ##history(conty)[[i]] <- history(cel) now through MIAME

    ival[, i] <- c(intensity(cel))
  }
  intensity(conty) <- ival
  #colnames(intensity(conty)) = filenames
  return(conty)
}

list.celfiles <-   function(...){
  files <- list.files(...)
  return(files[grep("\.[cC][eE][lL]\.gz$|\.[cC][eE][lL]$", files)])
}

###this is user friendly wrapper for read.affybatch
ReadAffy <- function(..., filenames=character(0),
                     widget=getOption("BioC")$affy$use.widgets,
                     compress=getOption("BioC")$affy$compress.cel,
                     celfile.path=getwd(),
                     sampleNames=NULL,
                     phenoData=NULL,
                     description=NULL,
                     notes="",
                     rm.mask=FALSE, rm.outliers=FALSE, rm.extra=FALSE,
                     verbose=FALSE) {
  
  ##first figure out filenames
  auxnames <- unlist(as.list(substitute(list(...)))[-1])

  if (widget){
    require(tkWidgets)
    widgetfiles <- fileBrowser(textToShow="Choose CEL files",
                               testFun=hasSuffix("[cC][eE][lL]"))
  }
  else
    widgetfiles <- character(0)
  
  filenames <- .Primitive("c")(filenames, auxnames, widgetfiles)
  
  if(length(filenames)==0) filenames <- list.celfiles(celfile.path,full.names=TRUE)
  
  if(length(filenames)==0) stop("No cel filennames specified and no cel files in specified directory:",celfile.path,"\n")
  
  
  ##now assign sampleNames if phenoData not given
  if(is.null(phenoData)){
    if(is.null(sampleNames)){
      if(widget){
        require(tkWidgets)
        tksn <- tkSampleNames(filenames=filenames)
        sampleNames <- tksn[,1]
        ##notice that a description of the files is ingored for now
        ##soon to go into MIAME
      }
      else{
        sampleNames <- sub("^/?([^/]*/)*", "", filenames, extended=TRUE)
      }
    }
    else{
      if(length(sampleNames)!=length(filenames)){
        warning("sampleNames not same length as filenames. Using filenames as sampleNames instead\n")
        sampleNames <- sub("^/?([^/]*/)*", "", filenames, extended=TRUE)
      }
    }
  }
  
  ##now get phenoData
  if(is.character(phenoData)) ##if character read file
    phenoData <- read.phenoData(filename=phenoData)
  else{
    if(class(phenoData)!="phenoData"){
      if(widget){
        require(tkWidgets)
        phenoData <- read.phenoData(sampleNames=sampleNames,widget=TRUE)
      }
      else
        phenoData <- read.phenoData(sampleNames=sampleNames,widget=FALSE)
    }
  }
  
  ##get MIAME information
  if(is.character(description)){
    description <- read.MIAME(filename=description,widget=FALSE)
  }
  else{
    if(class(description)!="MIAME"){
      if(widget){
        require(tkWidgets)
        description <- read.MIAME(widget=TRUE)
      }
      else
        description <- new("MIAME")
    }
  }
  
  ##MIAME stuff
  description@preprocessing$filenames <- filenames
  if(exists("tksn")) description@samples$description <- tksn[,2]
  description@preprocessing$affyversion <- library(help=affy)$info[[2]][[2]][2]

  ##and now we are ready to read cel files
  return(read.affybatch(filenames=filenames,
                        phenoData=phenoData,
                        description=description,
                        notes=notes,
                        compress=compress,
                        rm.mask=rm.mask,
                        rm.outliers=rm.outliers,
                        rm.extra=rm.extra,
                        verbose=verbose))
}








