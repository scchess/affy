

.First.lib <- function(libname, pkgname, where) {
  

  ## DEBUG flag
  assign("debug.affy123", T, envir=.GlobalEnv)
  
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
  
  require(methods, quiet=TRUE)
  require(modreg, quiet=TRUE)
  require(eda, quiet=TRUE)
  
  where <- match(paste("package:", pkgname, sep=""), search())

  if (debug.affy123) cat("-->detecting normalization methods from naming convention\n")
  ## this could move into the respective methods of Plob and Cel.container later
  assign("normalize.Cel.container.methods",
         substr(ls(where)[grep("normalize\.Cel\.container\.*", ls(where))], 25, 100),
         envir=as.environment(where))
  assign("normalize.Plob.methods",
         substr(ls(where)[grep("normalize\.Plob\.*", ls(where))], 16,100),
         envir=as.environment(where))

  if (debug.affy123) cat("-->detecting expression value methods from naming convention\n")
  ## the first one is deprecated (well... "should be"...)
  assign("generateExprSet.methods",
         substr(ls(where)[grep("generateExprVal\.method\.*", ls(where))], 24,100),
         envir=as.environment(where))
  assign("express.summary.stat.methods",
         substr(ls(where)[grep("generateExprVal\.method\.*", ls(where))], 24,100),
         envir=as.environment(where))

  if (debug.affy123) cat("-->detecting background correction methods from naming convention\n")
  ##assign("bg.correct.methods",
  ##       substr(ls(where)[grep("bg.correct\.*", ls(where))], 12,100),
  ##       envir=as.environment(where))
  assign("bg.correct.methods",
         ls(where)[grep("bg.correct\.*", ls(where))],
         envir=as.environment(where))
         
  
  ##if (debug.affy123) cat("-->initCdf\n")
  .initCdf(where)
  ##if (debug.affy123) cat("-->initCel\n")
  .initCel(where)
  ##if (debug.affy123) cat("-->initCel.container\n")
  .initCel.container(where)
  ##if (debug.affy123) cat("-->initPPSet\n")
  .initPPSet(where)
  ##if (debug.affy123) cat("-->initPPSet.container\n")
  .initPPSet.container(where)
  ##if (debug.affy123) cat("-->initPlob\n")
  .initPlob(where)


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





