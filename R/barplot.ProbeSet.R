barplot.ProbeSet <- function(height,
                             xlab="Probe pair",ylab="Intensity",
                             main=NA,
                             col.pm="red", col.mm="blue",
                             beside=TRUE, names.arg="pp",
                             ask = TRUE,
                             ...)
{

  opar <- par()$ask
  par(ask=ask)
  on.exit(par(ask=opar))

  if (names.arg == "pp") {
    names.arg <- seq(1, nrow(pm(height)))
  }

  col <- c(col.pm, col.mm)
  
  for (i in 1:ncol(pm(height))) {
    if (is.na(main))
      main <- paste(height@id, "(", i, ")")
    
    hh <- rbind(pm(height)[, i], mm(height)[, i])
      
    barplot(hh, xlab=xlab, ylab=ylab,
            main=main,
            col=col,
            beside=beside,
            names.arg=names.arg,
            ...)
  }
}
