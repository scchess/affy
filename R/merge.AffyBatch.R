merge.AffyBatch <- function(x, y, annotation=paste(x@annotation, y@annotation),
                            description=NULL,
                            notes=paste(x@notes, y@notes), ...) {

  adim <- dim(intensity(x))[1]

  if ((x@nrow != y@nrow) || (x@ncol != y@ncol))
    stop("cannot merge chips of different sizes !")

  if (x@cdfName != y@cdfName)
    warning("cdfName mismatch (using the cdfName of x)!")

  if (is.null(description))
    description <- paste("merged from ", deparse(x), " and ", deparse(y))
                         
  lx <- length(x)
  ly <- length(y)
  
  mlcel <- new("AffyBatch",
               intensity=array(NA, dim=c(adim, lx+ly)),
               cdfName=x@cdfName,
               nrow=x@nrow,
               ncol=x@ncol,
               nexp=lx+ly,
               ##sd=array(),
               chipNames=c(chipNames(x), chipNames(y)),
               ##outliers=vector("list", lx+ly),
               ##masks=vector("list", lx+ly),
               annotation=annotation,
               description=description,
               notes=notes,
               history=vector("list", lx+ly)
               )

  for (i in 1:lx) {
    intensity(mlcel)[, i] <- intensity(x)[, i]
    ##outliers(mlcel)[[i]] <- outliers(x)[[i]]
    ##masks(mlcel)[[i]] <- masks(x)[[i]]
    history(mlcel)[[i]] <- history(x)[[i]]
  }
  
  for (i in 1:ly) {
    intensity(mlcel)[, i+lx] <- intensity(y)[, i]
    ##outliers(mlcel)[[i+lx]] <- outliers(y)[[i]]
    ##masks(mlcel)[[i+lx]] <- masks(y)[[i]]
    history(mlcel)[[i+lx]] <- history(y)[[i]]
  }
  
  return(mlcel)
  
}
