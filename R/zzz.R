.initNormalize <- function(all.affy, env) {
  if (debug.affy123) cat("-->detecting normalization methods from naming convention\n")

  ## this could move into the respective methods of AffyBatch later

  assign("normalize.AffyBatch.methods",
         sub("^normalize\\.AffyBatch\\.", "",
             grep("^normalize.AffyBatch", all.affy, value = TRUE)),
         envir=env)
}

.initExpression <- function(all.affy, env) {
  if (debug.affy123) cat("-->detecting expression value methods from naming convention\n")

  ## the first one is deprecated (well... "should be"...)
  vals <- sub("^generateExprVal\\.method\\.", "",
              grep("^generateExprVal.method", all.affy, value = TRUE))
  assign("generateExprSet.methods", vals,
         envir=env)
  assign("express.summary.stat.methods", vals,
         envir=env)
}

.initBackgroundCorrect <- function(all.affy, env) {
  if (debug.affy123) cat("-->detecting background correction methods from naming convention\n")
  start <- nchar("bg.correct.")
  assign("bgcorrect.methods",
         sub("^bg\\.correct\\.", "",
             grep("^bg.correct", all.affy, value = TRUE)),
         envir=env)
       }

.initPmCorrect <- function(all.affy, env) {
  if (debug.affy123) cat("-->detecting pm correction methods from naming convention\n")
  
  assign("pmcorrect.methods",
         sub("^pmcorrect\\.", "",
             grep("^pmcorrect", all.affy, value = TRUE)),
         envir=env)
}

.setAffyOptions <- function(affy.opt=NA) {

  if (! any(is.na(affy.opt))) {
    if (class(affy.opt) != "BioCPkg")
      stop("obviously invalid package options !")

    BioC <- getOption("BioC")
    BioC$affy <- affy.opt
    options("BioC"=BioC)
    return()
  }

  ## add affy specific options
  ## (not unlike what is done in 'Biobase')
  if (is.null(getOption("BioC"))) {
    BioC <- list()
    class(BioC) <- "BioCOptions"
    options("BioC"=BioC)
  }

  probesloc.first <- list(what="environment", where=.GlobalEnv)
  probesloc.second <- list(what="libPath", where=NULL)
  probesloc.third <- list(what="data", where="affy")
  probesloc.fourth <- list(what="bioC", where=.libPaths()[1])


  ## default for the methods
  bgcorrect.method <- "mas"
  normalize.method <- "quantiles"
  pmcorrect.method <- "pmonly"
  summary.method <- "liwong"

  affy.opt <- list(compress.cdf=FALSE, compress.cel=FALSE,
                   use.widgets=FALSE,
                   probesloc = list(probesloc.first, probesloc.second,
                   probesloc.third, probesloc.fourth),
                   bgcorrect.method = bgcorrect.method,
                   normalize.method = normalize.method,
                   pmcorrect.method = pmcorrect.method,
                   summary.method = summary.method,
                   xy.offset = 0 ## this one is for temporary compatibility
                   ) 

  class(affy.opt) <- "BioCPkg"

  BioC <- getOption("BioC")
  BioC$affy <- affy.opt
  options("BioC"=BioC)
  ## ---
}

.onLoad <- function(libname, pkgname) {
   
#  where <- match(paste("package:", pkgname, sep=""), search())
  all.affy <- ls(environment(sys.function()))

   ##a place to store some variables that need to be accessed
   .affyInternalEnv <- new.env(parent=emptyenv())

   assign(".affyInternalEnv", .affyInternalEnv, 
     envir=topenv(parent.frame()))

  .initNormalize(all.affy, .affyInternalEnv)
  .initExpression(all.affy, .affyInternalEnv)
  .initBackgroundCorrect(all.affy, .affyInternalEnv)
  .initPmCorrect(all.affy, .affyInternalEnv)

  .setAffyOptions()

  if(.Platform$OS.type == "windows" && require(Biobase) && interactive()
        && .Platform$GUI ==  "Rgui"){
        addVigs2WinMenu("affy")
    }
}

