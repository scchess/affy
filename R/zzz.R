.First.lib <- function(libname, pkgname, where) {
  debug.affy123 <- T #DEBUG
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
  
  rm(debug.affy123)
  cacheMetaData(as.environment(where))
}

