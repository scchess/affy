.First.lib <- function(libname, pkgname, where) {
  
  ## DEBUG flag
  assign("debug.affy123", T, envir=.GlobalEnv)

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
  ## this will probably move into the respective methods of Plob and Cel.container
  ## note: inconsistency (my own fault) in the naming convention 'Cel' will probably
  ## have to move them to 'Cel.container'
  assign("normalize.Cel.container.methods", substr(ls(where)[grep("normalize\.Cel\.*", ls(where))],
                                                   15, 100),
                                                   envir=as.environment(where))
  assign("normalize.Plobs.methods", substr(ls(where)[grep("normalize\.Plob\.*", ls(where))],
                                          15,100),
                                          envir=as.environment(where))
  
  cacheMetaData(as.environment(where))
}






