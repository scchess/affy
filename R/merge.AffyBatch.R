merge.AffyBatch <- function(x, y, annotation=paste(x@annotation, y@annotation),
                            description=NULL,
                            notes=paste(x@notes, y@notes), ...) {

  adim <- dim(intensity(x))[1]

  if ((x@nrow != y@nrow) || (x@ncol != y@ncol))
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
             nrow=x@nrow,
             ncol=x@ncol,
             phenoData=phenodata,
             annotation=x@annotation,
             description=x@description, ##need to write a merge for MIAME
             notes=paste(x@notes,y@notes))
         )
}
