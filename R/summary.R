###these are summary functions they take matrices of probes x chips
###and return expression and se (when applicable)

##DEBUG: appending the se to the expression values in a same vector
##       is too much hackish (I think)... we need to think about something
##       better

avdiff <- function(x,constant=3){
  e <- apply(x,2,function(y){
    o <- order(y)
    yy <- y[-c(o[1],o[length(y)])] #take out biggest and smallest
    if(length(yy)<2)  # SK, some genes have only one probe
      mean(y)
    else
      mean(y[abs(y-mean(yy))<constant*sd(yy)])
  })
  c(e,rep(NA,length(e)))
  
}

li.wong <- function(data.matrix,remove.outliers=T,
                    normal.array.quantile=0.5,
                    normal.resid.quantile=0.9,
                    large.threshold=3,
                    large.variation=0.8,
                    outlier.fraction=0.14,
                    delta = 1e-06,maxit=50,outer.maxit=50,verbose=F){

  e <-  fit.li.wong(t(data.matrix),remove.outliers,normal.array.quantile,normal.resid.quantile,large.threshold,large.variation,outlier.fraction,delta,maxit,outer.maxit,verbose)
  c(e$theta,e$sigma.theta)
}


medianpolish <- function(x, ...){
  tmp <- medpolish(log2(x), trace.iter=F, ...)
  ##rough estimate
  sigma <- 1.483*median(abs(as.vector(tmp$residuals)))/sqrt(nrow(x))
  c(tmp$overall + tmp$col,rep(sigma, ncol(x)))
}


##DEBUG: to be moved to a proper place (eventually renamed)
generateExprVal.method.medianpolish <- function(matos, ...) {
  medianpolish(matos, ...)[1:ncol(matos)]
}

