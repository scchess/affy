merge.AffyBatch <- function(x, y, annotation=paste(annotation(x), annotation(y)),
                            description=NULL,
                            notes=paste(x@notes, y@notes), ...) {

  adim <- dim(intensity(x))[1]

  if ((nrow(x) != nrow(y)) || (ncol(x) != ncol(y)))
    stop("cannot merge chips of different sizes !")

  if (x@cdfName != y@cdfName)
    warning("cdfName mismatch (using the cdfName of x)!")

  if (is.null(description))
    description <- paste("merged")
                         
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
             description=description(x), ##need to write a merge for MIAME
             notes=paste(x@notes,y@notes))
         )
}
