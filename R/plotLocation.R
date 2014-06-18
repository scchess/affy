plotLocation <- function(x, col="green", pch=22, ...) {
  if (is.list(x)) {
    x <- cbind(unlist(lapply(x, function(x) x[,1])),
               unlist(lapply(x, function(x) x[,2])))
  }
  ## need to use nrow - x[,2] for correct y position.
  ## use image width to get nrow, which isn't ideal
  ## but follows assumption for this function that an
  ## image already exists
  nrow <- par("usr")[4]
  if(nrow == 1) stop(paste("\nYou must first generate an image of an array",
     "for this function to work!\n\n"), call. = FALSE)
  points(x[,1], nrow - x[,2]
         , pch=pch, col=col, ...)
}
