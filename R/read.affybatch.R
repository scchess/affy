#############################################################
##
## read.affybatch.R
##
## Adapted by B. M. Bolstad from read.affybatch in the affy
## package version 1.2.  The goal is a faster, less memory hungry
## ReadAffy. To do this we will shunt more work off to
## the c code.
##
## History
## Jun 13-15 Intial version
## Jun 16    Verbose flag passed to C routine
## Jun 17    New method for checking header of first cel
##           file.
## Jul 7     Added the function read.probematrix which
##           reads in PM, MM or both into matrices
## Sep 28    changed name from read.affybatch2 to read.affybatch
##           and cleaned up some old commented stuff
## Apr 13, 2004 - fixed problem in read.probematrix
##
#############################################################


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

  headdetails <- .Call("ReadHeader",filenames[[1]],compress, PACKAGE="affy")

  #print(headdetails)



  ##now we use the length
  dim.intensity <- headdetails[[2]]   ##dim(intensity(cel))
  ##and the cdfname as ref
  ref.cdfName <- headdetails[[1]]   #cel@cdfName

  if (verbose)
    cat(paste("instanciating an AffyBatch (intensity a ", prod(dim.intensity), "x", length(filenames), " matrix)...", sep=""))



  if (verbose)
    cat("done.\n")

  #### this is where the code changes from the original read.affybatch.
  #### what we will do here is read in from the 1st to the nth CEL file

  return(new("AffyBatch",
               exprs  = .Call("read_abatch",filenames,compress, rm.mask,
               rm.outliers, rm.extra, ref.cdfName,
               dim.intensity,verbose, PACKAGE="affy"),
               ##se.exprs = array(NaN, dim=dim.sd),
               cdfName    = ref.cdfName,   ##cel@cdfName,
               phenoData  = phenoData,
               nrow       = dim.intensity[1],
               ncol       = dim.intensity[2],
               annotation = cleancdfname(ref.cdfName, addcdf=FALSE),
               description= description,
               notes      = notes))
}





######################################################################################

read.probematrix <- function(..., filenames = character(0), phenoData = new("phenoData"),
    description = NULL, notes = "", compress = getOption("BioC")$affy$compress.cel,
    rm.mask = FALSE, rm.outliers = FALSE, rm.extra = FALSE, verbose = FALSE,which="pm"){

  auxnames <- as.list(substitute(list(...)))[-1]
  filenames <- .Primitive("c")(filenames, auxnames)

  match.arg(which,c("pm","mm","both"))

  if (verbose)
        cat(1, "reading", filenames[[1]], "to get header informatio")
    headdetails <- .Call("ReadHeader", filenames[[1]], compress, PACKAGE="affy")
    dim.intensity <- headdetails[[2]]
    ref.cdfName <- headdetails[[1]]

  Data <- new("AffyBatch", cdfName = ref.cdfName, annotation = cleancdfname(ref.cdfName,addcdf = FALSE))

  cdfInfo<- as.list(getCdfInfo(Data))
  cdfInfo <- cdfInfo[order(names(cdfInfo))]

  
  .Call("read_probeintensities", filenames,
        compress, rm.mask, rm.outliers, rm.extra, ref.cdfName,
        dim.intensity, verbose, cdfInfo,which, PACKAGE="affy")
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








