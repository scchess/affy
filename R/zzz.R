.First.lib <- function(libname, pkgname, where) {
  

  ## DEBUG flag
  assign("debug.affy123", T, envir=.GlobalEnv)
  if (debug.affy123) {
    cat(rep("*",13),"\n",sep="")
    cat("affy: development version (might be unstable)\n")
    cat(rep("*",13),"\n",sep="")
    cat("demo(affy.tour) will eventually work and give an overview...\n")
  }
  
  library.dynam("affy", pkgname, libname)
  require(methods)
  require(modreg)
  require(eda)
  where <- match(paste("package:", pkgname, sep=""), search())
  
  if (debug.affy123) cat("-->initCdf\n")
  .initCdf(where)
  if (debug.affy123) cat("-->initCel\n")
  .initCel(where)
  if (debug.affy123) cat("-->initCel.container\n")
  .initCel.container(where)
  if (debug.affy123) cat("-->initPPSet\n")
  .initPPSet(where)
  if (debug.affy123) cat("-->initPPSet.container\n")
  .initPPSet.container(where)
  if (debug.affy123) cat("-->initPlob\n")
  .initPlob(where)

  if (debug.affy123) cat("-->detecting normalization methods from naming convention\n")
  ## this could move into the respective methods of Plob and Cel.container later
  assign("normalize.Cel.container.methods", substr(ls(where)[grep("normalize\.Cel\.container\.*", ls(where))],
                                                   25, 100),
         envir=as.environment(where))
  assign("normalize.Plob.methods", substr(ls(where)[grep("normalize\.Plob\.*", ls(where))],
                                          16,100),
         envir=as.environment(where))

  if (debug.affy123) cat("-->detecting expression value methods from naming convention\n")
  assign("generateExprSet.methods", substr(ls(where)[grep("generateExprVal\.method\.*", ls(where))],
                                           23,100),
         
         envir=as.environment(where))
  if (debug.affy123) cat("-->detecting expression value methods from naming convention\n")
  assign("bg.correct.methods", substr(ls(where)[grep("bg.correct\.*", ls(where))],
                                           12,100),
         
         envir=as.environment(where))
  cacheMetaData(as.environment(where))

}






