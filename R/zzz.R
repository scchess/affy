.First.lib <- function(libname, pkgname, where) {
  library.dynam("affy", pkgname, libname)
  require(methods)
  where <- match(paste("package:", pkgname, sep=""), search())
  .initCel(where)
  .initCel.container(where)
  cacheMetaData(as.environment(where))
}

