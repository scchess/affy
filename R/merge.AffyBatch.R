merge.AffyBatch <- function(x, y, annotation=paste(annotation(x), annotation(y)),
                            description=NULL,
                            notes=paste(x@notes, y@notes), ...) {

  adim <- dim(intensity(x))[1]

  if ((nrow(x) != nrow(y)) || (ncol(x) != ncol(y)))
    stop("cannot merge chips of different sizes !")

  if (x@cdfName != y@cdfName)
    warning("cdfName mismatch (using the cdfName of x)!")

  if (is.null(description)){
    description <- new("MIAME")
    description@title <- "Created from merging two AffyBatches. No description was supplied. The description of the two original AffyBatches was not kept."
  }                       

  lx <- length(x)
  ly <- length(y)

  phenodata <- phenoData(x)
  pData(phenodata) <- rbind(pData(x),pData(y))
  return(new("AffyBatch",
             exprs=cbind(intensity(x),intensity(y)),
             cdfName=x@cdfName,
             nrow=nrow(x),
             ncol=ncol(x),
             phenoData=phenodata,
             annotation=x@annotation,
             description=description, ##need to write a merge for MIAME
             notes=paste("Merge from two AffyBatches with notes: 1)",
               x@notes,", and 2)",y@notes))
         )
}



