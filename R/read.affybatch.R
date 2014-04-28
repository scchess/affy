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
## Nov 15, 2005 - add functionality to read the
##                stddev values into the se.exprs slot (non-default behaviour)
##
## Jan 24, 2006 - JWM: added cdfname to allow for the use of non-standard mappings
## Mar 6, 2006 - change .Call to reference affyio. that is new location for parsing code
## Dec 12, 2006 - added checkCelFiles() to ensure all filenames are celfiles so unintended
##                arguments don't get passed in via ...
## Apr 19, 2013 - JWM: added warning and error messages for Gene ST and Exon ST arrays
## Sept 26, 2013 - naked .Call() to affyio replaced
##
#############################################################


read.affybatch <- function(..., filenames=character(0),
                           phenoData=new("AnnotatedDataFrame"),
                           description=NULL,
                           notes="",
                           compress = getOption("BioC")$affy$compress.cel,
                           rm.mask = FALSE, rm.outliers=FALSE, rm.extra=FALSE,
                           verbose = FALSE,sd=FALSE, cdfname = NULL) {

  auxnames <- unlist(list(...))
  filenames <- c(filenames, auxnames)
  checkValidFilenames(filenames)

  n <- length(filenames)
  pdata <- pData(phenoData)
  ## try to read sample names form phenoData. if not there use CEL
  ## filenames
  if(dim(pdata)[1] != n) {
    ## if empty pdata filename are samplenames
    warning("Incompatible phenoData object. Created a new one.\n")

    samplenames <- sub("^/?([^/]*/)*", "", filenames)
    pdata <- data.frame(sample=1:n, row.names=samplenames)
    phenoData <- new("AnnotatedDataFrame",
                     data=pdata,
                     varMetadata=data.frame(
                       labelDescription="arbitrary numbering",
                       row.names="sample"))
  }
  else samplenames <- rownames(pdata)

  if (is.null(description))
    {
      description <- new("MIAME")
      preproc(description)$filenames <- filenames
      preproc(description)$affyversion <- library(help=affy)$info[[2]][[2]][2]
    }
  if (length(notes)==0) notes(description) <- notes
  ## read the first file to see what we have
  if (verbose) cat(1, "reading",filenames[[1]],"...")

  headdetails <- read.celfile.header(as.character(filenames[[1]]))
  ##now we use the length
  dim.intensity <- headdetails[[2]]   ##dim(intensity(cel))
  ##and the cdfname as ref
  ref.cdfName <- headdetails[[1]]   #cel@cdfName

  if(length(grep("gene1[01]st", cleancdfname(ref.cdfName))) == 1)
      warning(paste0("\n\nThe affy package can process data from the Gene ST 1.x series of arrays,\n",
                    "but you should consider using either the oligo or xps packages, which are specifically\n",
                     "designed for these arrays.\n\n"), call. = FALSE)
  if(length(grep("gene2[01]st|ex[1-2][0-1]st|hta20", cleancdfname(ref.cdfName))) == 1)
      stop(paste0("\n\nThe affy package is not designed for this array type.\n",
                   "Please use either the oligo or xps package.\n\n"), call. = FALSE)
              

  scandates <-
    sapply(seq_len(length(filenames)), function(i) {
             sdate <- read.celfile.header(filenames[i], info = "full")[["ScanDate"]]
             if (is.null(sdate) ||length(sdate) == 0 ) NA_character_ else sdate
           })
  protocol <-
    new("AnnotatedDataFrame",
        data=data.frame("ScanDate"=scandates, row.names=sampleNames(phenoData),
                        stringsAsFactors=FALSE),
        dimLabels=c("sampleNames", "sampleColumns"))

  ## allow for non-standard cdfs
  if(is.null(cdfname))
    cdfname <- ref.cdfName

  if (verbose)
    cat(paste("instantiating an AffyBatch (intensity a ", prod(dim.intensity), "x", length(filenames), " matrix)...", sep=""))

  if (verbose)
    cat("done.\n")

  ## Change sampleNames to be consistent with row.names of phenoData
  ## object

  exprs <-  affyio::read_abatch(filenames, rm.mask,
               rm.outliers, rm.extra, ref.cdfName,
               dim.intensity[c(1,2)],verbose)
  colnames(exprs) <- samplenames
  
  #### this is where the code changes from the original read.affybatch.
  #### what we will do here is read in from the 1st to the nth CEL file
  if (!sd){
    return(new("AffyBatch",
               exprs  = exprs,
               ##se.exprs = array(NaN, dim=dim.sd),
               cdfName    = cdfname,   ##cel@cdfName,
               phenoData  = phenoData,
               nrow       = dim.intensity[2],##["Rows"],
               ncol       = dim.intensity[1],##["Cols"],
               annotation = cleancdfname(cdfname, addcdf=FALSE),
               protocolData  = protocol,
               description= description,
               notes      = notes))
  } else {
    return(new("AffyBatch",
               exprs  = exprs,
               se.exprs =  affyio::read_abatch_stddev(filenames, rm.mask,
               rm.outliers, rm.extra, ref.cdfName,
               dim.intensity[c(1,2)],verbose),
               cdfName    = cdfname,   ##cel@cdfName,
               phenoData  = phenoData,
               nrow       = dim.intensity[2],##["Rows"],
               ncol       = dim.intensity[1],##["Cols"],
               annotation = cleancdfname(cdfname, addcdf=FALSE),
               protocolData  = protocol,
               description= description,
               notes      = notes))
  }
}





######################################################################################

read.probematrix <- function(..., filenames = character(0), phenoData = new("AnnotatedDataFrame"),
                             description = NULL, notes = "", compress = getOption("BioC")$affy$compress.cel,
                             rm.mask = FALSE, rm.outliers = FALSE, rm.extra = FALSE, verbose = FALSE,which="pm",
                             cdfname = NULL){

  auxnames <- unlist(list(...))
  filenames <- c(filenames, auxnames)

  which <- match.arg(which,c("pm","mm","both"))

  if (verbose)
    cat(1, "reading", filenames[[1]], "to get header information\n")
  headdetails <- read.celfile.header(as.character(filenames[[1]]))
  ref.cdfName <- headdetails[[1]]
  cleaned.cdfName <- cleancdfname(ref.cdfName, addcdf = FALSE)
  ## Allow for usage of alternative cdfs
  if(is.null(cdfname))
    Data <- new("AffyBatch", cdfName = ref.cdfName, annotation = cleaned.cdfName)
  else
    Data <- new("AffyBatch", cdfName = cdfname, annotation = cleaned.cdfName)
  
  cdfInfo <- as.list(getCdfInfo(Data))
  cdfInfo <- cdfInfo[order(names(cdfInfo))]

  read.celfile.probeintensity.matrices(filenames = filenames,
                                       cdfInfo = cdfInfo,
                                       rm.mask = rm.mask,
                                       rm.outliers = rm.outliers,
                                       rm.extra = rm.extra,
                                       verbose = verbose,
                                       which = which)
}


list.celfiles <-   function(...){
  files <- list.files(...)
  return(files[grep("\\.[cC][eE][lL]\\.gz$|\\.[cC][eE][lL]$", files)])
}

AllButCelsForReadAffy <- function(..., filenames=character(0),
                                  widget=getOption("BioC")$affy$use.widgets,
                                  celfile.path=NULL,
                                  sampleNames=NULL,
                                  phenoData=NULL,
                                  description=NULL){

  ##first figure out filenames
  auxnames <- unlist(as.list(substitute(list(...)))[-1])

  if (widget){
    require(tkWidgets)
    widgetfiles <- fileBrowser(textToShow="Choose CEL files",
                               testFun=hasSuffix("[cC][eE][lL]|[cC][eE][lL].gz"))
  }
  else{
    widgetfiles <- character(0)
  }

  if(!is.null(celfile.path)){
    auxnames <- file.path(celfile.path, auxnames)
    filenames <- file.path(celfile.path, filenames)
  }

  filenames <- c(filenames, auxnames, widgetfiles)

  if(length(filenames)==0){
    if(is.null(celfile.path)) celfile.path <- getwd()
    filenames <- list.celfiles(celfile.path,full.names=TRUE)
  }
  if(length(filenames)==0) stop("No cel filennames specified and no cel files in specified directory:",celfile.path,"\n")

  if(is.null(sampleNames)){
    sampleNames <- sub("^/?([^/]*/)*", "", filenames)
  }
  else{
    if(length(sampleNames)!=length(filenames)){
      warning("sampleNames not same length as filenames. Using filenames as sampleNames instead\n")
      sampleNames <- sub("^/?([^/]*/)*", "", filenames)
    }
  }

  chkSn <- function(filenames, samplenames){
      fntest <- sub("^/?([^/]*/)*", "", filenames)
      if(all(fntest %in% samplenames)){
          filenames <<- filenames[match(samplenames, fntest)]
      } else {
          warning(paste("The cel file names are not conformable to the sample names",
                        "provided in the phenoData object.\nNote that the sample names",
                        "take precedence, so will be used in the sampleNames slot of the",
                        "resulting ExpressionSet.\nPlease ensure that the ordering of the",
                        "phenoData object matches the order of the file names",
                        "(use list.celfiles() to test), as no re-ordering is possible.\n\n"),
                  call. = FALSE)
      }
  }
  
  if(is.character(phenoData)) { 
    ## if character, read file 
    if(length(phenoData)!=1) stop(sprintf("'phenoData' must be of length 1, but is %d.", length(phenoData)))
    phenoData <- read.AnnotatedDataFrame(filename=phenoData)
    sampleNames <- sampleNames(phenoData)
    chkSn(filenames, sampleNames)
  } else if(is.data.frame(phenoData)) {
    ## if data.frame, coerce
    phenoData <-  as(phenoData, "AnnotatedDataFrame")
    sampleNames <- sampleNames(phenoData)
    chkSn(filenames, sampleNames)
  } else if(is.null(phenoData)) {
    phenoData <- new("AnnotatedDataFrame",
                     data = data.frame(sample=seq_along(sampleNames),
                                       row.names=sampleNames),
                     varMetadata = data.frame(labelDescription="arbitrary numbering",
                                    row.names=names(pData)))
  } else if (!is(phenoData, "AnnotatedDataFrame")) {
    stop(sprintf("'phenoData' must be of class 'AnnotatedDataFrame', but is %s.", class(phenoData)))
  }
  
  ##get MIAME information
  if(is.character(description)){
    description <- read.MIAME(filename=description,widget=FALSE)
  }
  else{
      if (! is(description, "MIAME")) {
          description <- new("MIAME")
      }
  }

  ##MIAME stuff
  description@preprocessing$filenames <- filenames
  description@preprocessing$affyversion <- library(help=affy)$info[[2]][[2]][2]

  return(list(filenames   = filenames,
              phenoData   = phenoData,
              sampleNames = sampleNames,
              description = description))
}

###this is user friendly wrapper for read.affybatch
ReadAffy <- function(..., filenames=character(0),
                     widget=getOption("BioC")$affy$use.widgets,
                     compress=getOption("BioC")$affy$compress.cel,
                     celfile.path=NULL,
                     sampleNames=NULL,
                     phenoData=NULL,
                     description=NULL,
                     notes="",
                     rm.mask=FALSE, rm.outliers=FALSE, rm.extra=FALSE,
                     verbose=FALSE,sd=FALSE, cdfname = NULL) {

  l <- AllButCelsForReadAffy(..., filenames=filenames,
                             widget=widget,
                             celfile.path=celfile.path,
                             sampleNames=sampleNames,
                             phenoData=phenoData,
                             description=description)

  ##and now we are ready to read cel files
  ret <- read.affybatch(filenames=l$filenames,
                        phenoData=l$phenoData,
                        description=l$description,
                        notes=notes,
                        compress=compress,
                        rm.mask=rm.mask,
                        rm.outliers=rm.outliers,
                        rm.extra=rm.extra,
                        verbose=verbose,sd=sd,cdfname=cdfname)

  sampleNames(ret) <- l$sampleNames
  return(ret)
}

checkValidFilenames <- function(filenames) {
    ## Returns TRUE if filenames is a character vector containing
    ## paths to files that exist (directories don't count).
    ## A suitable error message is printed via stop() if invalid
    ## file names are encountered.
    if (!is.character(filenames))
      stop(strwrap(paste("file names must be specified using a character",
                         "vector, not a", sQuote(typeof(filenames)))),
           call.=FALSE)

    if (length(filenames) == 0)
      stop("no file names provided")

    if (any(sapply(filenames, nchar) < 1))
      stop("empty file names are not allowed")

    finfo <- file.info(filenames)
    whBad <- sapply(finfo[["isdir"]], function(x) !identical(FALSE, x))
    if (any(whBad)) {
        msg <- paste("the following are not valid files:\n",
                     paste("  ", filenames[whBad], collapse="\n"))
        stop(msg, call.=FALSE)
    }
    TRUE
}
