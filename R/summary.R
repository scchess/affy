###these are summary functions they take matrices of probes x chips
###and return expression and se (when applicable)
avdiff <- function(x,constant=3){
  e <- apply(x,2,function(y){
    o <- order(y)
    yy <- y[-c(o[1],o[length(y)])] #take out biggest and smallest
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


medianpolish <- function(x,...){
  tmp <- medpolish(log2(x),trace.iter=F,...)
  ##rough estimate
  sigma_1.483*median(abs(as.vector(tmp$residuals)))/sqrt(nrow(x))
  c(tmp$overall + tmp$col,rep(sigma,ncol(x)))
}




