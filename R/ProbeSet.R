.initProbeSet <- function(where) {
 ## A ProbeSet holds probe values for a probe pairs set(*) accross a batch of experiments.
  ## methods 'express.summary.stat' returns of expression value per experiement in the
  ## batch, and 'bg.correct' does background correction (in some sense... the MM probes
  ## were created to measure unspecific hybridization. People thought that doing
  ## PM - MM would remove background noise. The method 'bg.correct' accepts extra parameters
  ## through '...' (can be used to pass background correction parameters common to different
  ## ProbeSet)
  ##
  ## - 
  ## (*) : a probe pair set is the set of probes pairs(**) related to an affyid. Generally a
  ##       a probe pair set has 20 elements.
  ## (**): a probe pair (or atom) is a pair of PM/MM values
  ##
  
  if (debug.affy123) cat("-->initProbeSet\n")
  
  setClass("ProbeSet",
           representation(id="character", pm="matrix", mm="matrix"),
           prototype=list(pm=matrix(), mm=matrix()),
           where=where)
  
  ###we need a show
  ## --> not only... barplot and plot for PPSet should have been ported too...

  setMethod("show", "ProbeSet",
            function(object) {
              cat("ProbeSet object:\n")
              cat("  id=", object@id, "\n", sep="")
              cat("  pm=", nrow(object@pm), "probes\n")
            }, where=where)
  
  ##DEBUG: what to do with that ?
  ## --> with what ?
  
  if( !isGeneric("colnames") )
    #setGeneric("colnames", function(x, do.NULL, prefix)
    setGeneric("colnames", where=where)
  
  ##for consistency also use sampleNames
  if( !isGeneric("sampleNames") )
    setGeneric("sampleNames", function(object)
               standardGeneric("sampleNames"), where=where)
  setMethod("sampleNames", "ProbeSet",
            function(object) colnames(object), where=where)
            
  setMethod("colnames", signature(x="ProbeSet"),
            function(x ,do.NULL=FALSE, prefix="row") {
              
              cnames<-colnames(pm(x))
              
              if (is.null(cnames)) {
                
                if (do.NULL) {
                  warning("No column names for ProbeSet")
                }
                else {
                  cnames <- paste(prefix, 1:ncols(x@pm))
                }
                
              } 
              return(cnames)
            },
            where=where)
  
  ## pm
  if( !isGeneric("pm") )
    setGeneric("pm", function(object) standardGeneric("pm"), where=where)
  
  setMethod("pm", "ProbeSet", function(object) object@pm)
  
  if( !isGeneric("pm<-") )
   setGeneric("pm<-", function(object, value) standardGeneric("pm<-"), where=where)
  
  setReplaceMethod("pm", signature=c("ProbeSet", "matrix"),
                   function(object, value) {
                     if (sum(dim(value) != dim(object@mm)) != 2)
                       stop("dimension mismatch between 'pm' and 'mm'")
                     object@pm <- value
                   })

  ## mm
  if( !isGeneric("mm") )
    setGeneric("mm", function(object) standardGeneric("mm"), where=where)
  setMethod("mm", "ProbeSet", function(object) object@mm)
  
  if( !isGeneric("mm<-") )
   setGeneric("mm<-", function(object, value) standardGeneric("mm<-"), where=where)
  
  setReplaceMethod("mm", signature=c("ProbeSet", "matrix"),
                   function(object, value) {
                     if (sum(dim(value) == dim(object@mm)) != 2)
                       stop("dimension mismatch between 'pm' and 'mm'")
                     object@mm <- value
                   })

  ## method express.summary.stat
  if( !isGeneric("express.summary.stat"))
    setGeneric("express.summary.stat", function(x, pmcorrect, method, ...)
               standardGeneric("express.summary.stat"), where=where)
  
  setMethod("express.summary.stat",signature(x="ProbeSet",  pmcorrect="character", method="character"),
            function(x, method, pmcorrect, pmcorrect.param=list(),
                     param.method=list()) {
            
              ## simple for system to let one add background correction methods
              ## relies on naming convention
              methodname <- paste("generateExprVal.method.", method, sep="")
              
              if (! exists(methodname))
                stop(paste("Unknown method (cannot find function", methodname, ")"))
              
              ## NOTE: this could change...
              #m <- do.call(bg.correct, c(alist(x@pm, x@mm), param.bg.correct))
              pm.corrected <- do.call(pmcorrect, c(alist(x@pm, x@mm), param.bg.correct))
              r <- do.call(methodname, c(alist(pm.corrected), param.method))
              
              ##DEBUG: name stuff to sort
              #names(r) <- names(allprobes)
              
              return(r)
            },
            where=where)
  ##barplot... to pass check
   if (debug.affy123) cat("--->probeset barplot\n")

  if( !isGeneric("barplot") )
    setGeneric("barplot",where=where)
  setMethod("barplot",signature(height="ProbeSet"),function(height,...) barplot.ProbeSet(height,...),where=where)
 
} 

###this is now in AffyBatch.R cause background must be performed on array
### method bg.correct 
#  if( !isGeneric("bg.correct") )
#    setGeneric("bg.correct", function(x, method, ...)
#               standardGeneric("bg.correct"), where=where)
  
##we no longer use bg correction on probe set. rather we do it on whole array
#   setMethod("bg.correct", signature(x="ProbeSet", method="character"),
#             function(x, method, ...) {

#               ## simple for system to let one add background correction methods
#               ## relies on naming convention
              
#               methodname <- paste("bg.correct.", method, sep="")
              
#               if (! exists(methodname))
#                 stop(paste("Unknown method (cannot find function", methodname, ")"))
              
#               r <- do.call(methodname, alist(x@pm, x@mm, ...))
              
#               return(r)
#             }, where=where)
# }

######no need for the commented out below since its no longer container
#   ## extract/set a PPSet in the container
#   setMethod("[[", "ProbeSet",
#             function(x, i, j) {
#               new("ProbeSet", probes=data.frame(pm=x@pm[, i], mm=x@mm[, i]),
#                   name=as.character(dimnames(x@pm)[[2]][i]))
#             },
#             where=where)
  
#   setReplaceMethod("[[", "ProbeSet",
#                    function(x, i, j, ..., value) {
#                      x@pm[, i] <- probes(value)$pm
#                      x@mm[, i] <- probes(value)$mm
#                      if (is.null(dimnames(x@pm))) {
#                        dimnames(x@pm) <- vector("list", 2)
#                        dimnames(x@pm)[[1]] <- vector("character", nrow(x@pm))
#                        dimnames(x@pm)[[2]] <- vector("character", ncol(x@pm))
#                      }
#                      dimnames(x@pm)[[2]][i] <- probeName(value)
#                      return(x)
#                    },
#                    where=where)
  

## while there is no generic for this one...
#names.ProbeSet <- function(x) {
#  return(as.character(lapply(x@x, function(y) y@name)))
#}


