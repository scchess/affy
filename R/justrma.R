## Sept 11, 2003 - justRMA calls just.rma2
### A user friendly wrapper for just.rma
justRMA <- function(..., filenames=character(0),
                     widget=getOption("BioC")$affy$use.widgets,
                     compress=getOption("BioC")$affy$compress.cel,
                     celfile.path=getwd(),
                     sampleNames=NULL,
                     phenoData=NULL,
                     description=NULL,
                     notes="",
                     rm.mask=FALSE, rm.outliers=FALSE, rm.extra=FALSE,
                     hdf5=FALSE, hdf5FilePath=NULL,verbose=FALSE,
                     normalize=TRUE, background=TRUE,
                     bgversion=2, destructive=FALSE){
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

  if(length(filenames)==0) stop("No cel filenames specified and no cel files in specified directory:",celfile.path,"\n")


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
  return(just.rma(filenames=filenames,
                        phenoData=phenoData,
                        description=description,
                        notes=notes,
                        compress=compress,
                        rm.mask=rm.mask,
                        rm.outliers=rm.outliers,
                        rm.extra=rm.extra,
                        #hdf5=hdf5, ## took these two out b/c I am not sure if hdf5 should be used
                        #hdf5FilePath=hdf5FilePath,
                        verbose=verbose,
                        normalize=normalize,
                        background=background,
                        bgversion=bgversion,
				destructive=destructive))
}




###########################################################################################
#
# this function uses a different parsing routine
# It was added Jul 7, 2003 by B. M. Bolstad
#
###########################################################################################

just.rma <- function(..., filenames=character(0),
                     phenoData=new("phenoData"),
                     description=NULL,
                     notes="",
                     compress=getOption("BioC")$affy$compress.cel,
                     rm.mask=FALSE, rm.outliers=FALSE, rm.extra=FALSE,
                     verbose=FALSE, background=TRUE, normalize=TRUE,
                     bgversion=2, destructive=FALSE) {

  auxnames <- as.list(substitute(list(...)))[-1]
  filenames <- .Primitive("c")(filenames, auxnames)

  n <- length(filenames)

  ## error if no file name !
  if (n == 0)
    stop("No file name given !")

  pdata <- pData(phenoData)
  ##try to read sample names form phenoData. if not there use CEL filenames
  if(dim(pdata)[1]!=n){#if empty pdata filename are samplenames
    warning("Incompatible phenoData object. Created a new one.\n")

    samplenames <- gsub("^/?([^/]*/)*", "", unlist(filenames), extended=TRUE	)
    pdata <- data.frame(sample=1:n,row.names=samplenames)
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
  ##if (verbose) cat(1, "reading",filenames[[1]],"...")

  ## get information from cdf environment

  headdetails <- .Call("ReadHeader", filenames[[1]], compress, PACKAGE="affy")
  dim.intensity <- headdetails[[2]]
  cdfname <- headdetails[[1]]

  tmp <- new("AffyBatch",
             cdfName=cdfname,
             annotation=cleancdfname(cdfname, addcdf=FALSE))
  pmIndex <- pmindex(tmp)
  probenames <- rep(names(pmIndex), unlist(lapply(pmIndex,length)))
  pmIndex <- unlist(pmIndex)

  ## read pm data into matrix

  probeintensities <- read.probematrix(filenames=filenames)

  ##pass matrix of pm values to rma

  ngenes <- length(geneNames(tmp))

  ##background correction
  bg.dens <- function(x){density(x,kernel="epanechnikov",n=2^14)}

  #if(destructive){
    exprs <-
    .Call("rma_c_complete",probeintensities$pm,probeintensities$pm,probenames,ngenes,body(bg.dens),new.env(),normalize,background,bgversion,
    PACKAGE="affy")
  #}else{
  #  exprs <- .Call("rma_c_complete_copy",probeintensities$pm,probeintensities$pm,probenames,ngenes,body(bg.dens),new.env(),normalize,background,bgversion)
  #}
  colnames(exprs) <- samplenames
  se.exprs <- array(NA, dim(exprs))

  annotation <- annotation(tmp)

  new("exprSet", exprs = exprs, se.exprs = se.exprs, phenoData = phenoData,
      annotation = annotation, description = description, notes = notes)
}





