.initNormalize <- function(where, all.affy) {
  if (debug.affy123) cat("-->detecting normalization methods from naming convention\n")

  #all.affy <- ls(where)
  
  ## this could move into the respective methods of Plob and Cel.container later
  assign("normalize.Cel.container.methods",
         substr(all.affy[grep("normalize\.Cel\.container\.*", all.affy)], 25, 100),
         envir=as.environment(where))
  assign("normalize.Plob.methods",
         substr(all.affy[grep("normalize\.Plob\.*", all.affy)], 16,100),
         envir=as.environment(where))
  assign("normalize.AffyBatch.methods",
         substr(all.affy[grep("normalize\.AffyBatch\.*", all.affy)], 21, 100),
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
  assign("bg.correct.methods",
         all.affy[grep("bg.correct\.*", all.affy)],
         envir=as.environment(where))
       }

.First.lib <- function(libname, pkgname, where) {
  

  where <- match(paste("package:", pkgname, sep=""), search())
  all.affy <- ls(where)
  
  ## DEBUG flag
  assign("debug.affy123", FALSE, envir=as.environment(where))

  ## R CMD check surprise...
  ## but the following trick did cheat R check... :(
  ##assign("T", TRUE, envir=as.environment(where))
  ##assign("F", FALSE, envir=as.environment(where))
  
  if (debug.affy123) {
    cat(rep("*",13),"\n",sep="")
    cat("affy: development version (...this is usable...)\n")
    cat("read the file NEWS to know what is changing.\n")
    cat(rep("*",13),"\n",sep="")
    cat("demo(affy.tour) will eventually work and give an overview...\n")
    cat(rep("*",13),"\n",sep="")
    cat("If you have tkWidgets and rhdf5 installed, try\nexpresso(widget=T, hdf5=T)\n for a real thrill...\n")
  }
  
  library.dynam("affy", pkgname, libname)
  
  require(methods, quietly=TRUE)
  require(modreg, quietly=TRUE)
  require(eda, quietly=TRUE)
  
  
  .initNormalize(where, all.affy)
  .initExpression(where, all.affy)
  .initBackgroundCorrect(where, all.affy)
  
  .initCdf(where)
  .initCel(where)
  .initPPSet(where)
  .initPPSet.container(where)
  .initAffyBatch(where)
  .initPlob(where)
  .initCel.container(where)

  ## add affy specific options
  ## (not unlike what is done in 'Biobase')
  if (is.null(getOption("BioC"))) {
    BioC <- list()
    class(BioC) <- "BioCOptions"
    options("BioC"=BioC)
  }
  
  ##affy$urls <- list( bioc = "http://www.bioconductor.org")
  affy <- list(compress.cdf=FALSE, compress.cel=FALSE,
               probesloc.what="package", probesloc.where=NULL, probesloc.autoload=TRUE)
  class(affy) <- "BioCPkg"
  
  BioC <- getOption("BioC")
  BioC$affy <- affy
  options("BioC"=BioC)
  ## ---
  
  cacheMetaData(as.environment(where))

}

.Last.lib <- function(libpath) {
  options(BioC)$affy <- NULL
  dyn.unload(file.path(libpath, "libs",
                       paste("affy", .Platform$"dynlib.ext", sep="")))
  .Dyn.libs <- .Dyn.libs[- which(.Dyn.libs == "affy")]
}





