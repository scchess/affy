######################################################
#
# rma - RMA interface to c code
#
# the RMA method implemented in c code
#
# this code serves as interface to the c code.
# currently
# implemented (version 0.25) background correction
#
# Background correction code has been added.
#
# note this function does not leave the supplied
# AffyBatch unchanged if you select DESTRUCTIVE=TRUE. this is 
# for memory purposes but can be quite
# dangerous if you are not careful. Use destructive=FALSE if this is
# deemed likely to be a problem.
#
########################################################

rma <- function(object,subset=NULL, verbose=TRUE, destructive = FALSE,...){

  rows <- length(probeNames(object))
  cols <- length(object)
 
  ngenes <- length(geneNames(object))
  
  #background correction
  bg.dens <- function(x){density(x,kernel="epanechnikov",n=2^14)}

  if (destructive){
  	exprs <- .Call("rma_c_complete",pm(object),mm(object),probeNames(object),ngenes,body(bg.dens),new.env())
  } else {
	exprs <- .Call("rma_c_complete_copy",pm(object),mm(object),probeNames(object),ngenes,body(bg.dens),new.env())
  }
  colnames(exprs) <- sampleNames(object)
  se.exprs <- array(NA, dim(exprs)) # to be fixed later, besides which don't believe much in nominal se's with medianpolish
  
  phenodata <- phenoData(object)
  annotation <- annotation(object)
  description <- description(object) 
  notes <- notes(object)
  
  new("exprSet", exprs = exprs, se.exprs = se.exprs, phenoData = phenodata, 
       annotation = annotation, description = description, notes = notes)
}
