.initNormalize <- function(where, all.affy) {
  if (debug.affy123) cat("-->detecting normalization methods from naming convention\n")
  
  ## this could move into the respective methods of AffyBatch later

  start <- nchar("normalize.AffyBatch.")
  assign("normalize.AffyBatch.methods",
         substr(all.affy[grep("normalize\.AffyBatch\.*", all.affy)], start+1, 100),
         envir=as.environment(where)) 
}

.initExpression <- function(where, all.affy) {
  if (debug.affy123) cat("-->detecting expression value methods from naming convention\n")
  
  ## the first one is deprecated (well... "should be"...)
  assign("generateExprSet.methods",
         substr(all.affy[grep("generateExprVal\.method\.*", all.affy)], 24,100),
         envir=as.environment(where))
  assign("express.summary.stat.methods",
         substr(all.affy[grep("generateExprVal\.method\.*", all.affy)], 24,100),
         envir=as.environment(where))
}

.initBackgroundCorrect <- function(where, all.affy) {
  if (debug.affy123) cat("-->detecting background correction methods from naming convention\n")
  ##assign("bg.correct.methods",
  ##       substr(ls(where)[grep("bg.correct\.*", ls(where))], 12,100),
  ##       envir=as.environment(where))
  start <- nchar("bg.correct.")
  assign("bgcorrect.methods",
         substr(all.affy[grep("bg\.correct\.*", all.affy)], start+1, 100),
         envir=as.environment(where))
       }

.initPmCorrect <- function(where, all.affy) {
  if (debug.affy123) cat("-->detecting pm correction methods from naming convention\n")
  start <- nchar("pmcorrect.")
  assign("pmcorrect.methods",
         substr(all.affy[grep("pmcorrect\.*", all.affy)], start+1, 100),
         envir=as.environment(where))
}

.initMapCdfName <- function(where) {
  filepath <- file.path(.path.package("affy"), "data", "mapCdfName.tab")
  mapCdfName <- read.table(filepath, colClasses=rep("character", 3), quote="\"", sep="\t", comment="#", row.names=NULL, header=TRUE)
  assign("mapCdfName", mapCdfName, envir=as.environment(where))
}

.First.lib <- function(libname, pkgname, where) {
  

  where <- match(paste("package:", pkgname, sep=""), search())
  all.affy <- ls(where)
  
  ## DEBUG flag
  assign("debug.affy123", TRUE, envir=as.environment(where))
  ##assign("debug.affy123", FALSE, envir=as.environment(where))
  
  message <- TRUE
  
  if (message) {
    cat(rep("*",13),"\n",sep="")
    cat("affy: development version\n")
    cat(rep("*",13),"\n",sep="")
    cat("The package is under major changes.\n")
    cat("unpack the package and read the file NEWS to know more....\n")
    cat("The draft for the new vignette (called 'affy2') is distributed with the pacakge\n")
    cat(rep("*",13),"\n",sep="")
    cat("demo(affy.tour) will eventually work and give an overview...\n")
    cat(rep("*",13),"\n",sep="")
  }
  
  library.dynam("affy", pkgname, libname)
  
  require(Biobase, quietly=TRUE) ##Biobase uses methods
  require(modreg, quietly=TRUE)
  require(eda, quietly=TRUE)

  ##i was having troulbes, and changing where to
  ###match(paste("package:", pkgname, sep=""), search()) fixed.. thanx to RG
  
  .initNormalize(match(paste("package:", pkgname, sep=""), search()), all.affy)
  .initExpression(match(paste("package:", pkgname, sep=""), search()), all.affy)
  .initBackgroundCorrect(match(paste("package:", pkgname, sep=""), search()), all.affy)
  .initPmCorrect(match(paste("package:", pkgname, sep=""), search()), all.affy)
  .initMapCdfName(match(paste("package:", pkgname, sep=""), search()))
  .initCdf(match(paste("package:", pkgname, sep=""), search()))
  .initCel(match(paste("package:", pkgname, sep=""), search()))
  .initAffyBatch(match(paste("package:", pkgname, sep=""), search()))
  .initProbeSet(match(paste("package:", pkgname, sep=""), search()))


  ## add affy specific options
  ## (not unlike what is done in 'Biobase')
  if (is.null(getOption("BioC"))) {
    BioC <- list()
    class(BioC) <- "BioCOptions"
    options("BioC"=BioC)
  }
  
  ##affy$urls <- list( bioc = "http://www.bioconductor.org")

  probesloc.first <- list(what="environment", where=.GlobalEnv)
  probesloc.second <- list(what="package", where=NULL, probesloc.autoload=TRUE)
  probesloc.third <- list(what="data", where="affy")

  
  ## i added use.widgets=FALSE. Shuold it be true?
  ## --> I do not think so. Let's keep it FALSE. 
  affy <- list(compress.cdf=FALSE, compress.cel=FALSE, use.widgets=FALSE,
               probesloc = list(probesloc.first, probesloc.second, probesloc.third))
  class(affy) <- "BioCPkg"
  
  BioC <- getOption("BioC")
  BioC$affy <- affy
  options("BioC"=BioC)
  ## ---

  cacheMetaData(as.environment(where))

}

.Last.lib <- function(libpath) {
  options("BioC")$affy <- NULL
  dyn.unload(file.path(libpath, "libs",
                       paste("affy", .Platform$"dynlib.ext", sep="")))
  .Dyn.libs <- .Dyn.libs[- which(.Dyn.libs == "affy")]
}





