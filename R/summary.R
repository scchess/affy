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

li.wong <- function(x,remove.outliers=T,
                    normal.array.quantile=0.5,
                    normal.resid.quantile=0.9,
                    large.threshold=3,
                    large.variation=0.8,
                    outlier.fraction=0.14,
                    delta = 1e-06,maxit=50,outer.maxit=50,verbose=F){

  e <-  fit.li.wong(t(x),remove.outliers,normal.array.quantile,normal.resid.quantile,large.threshold,large.variation,outlier.fraction,delta,maxit,outer.maxit,verbose)
  c(e$theta,e$sigma.theta)
}

medianpolish <- function(x,...){
  tmp <- medpolish(log2(x),trace.iter=F,...)
  c(tmp$overall + tmp$col,rep(NA,ncol(x)))
}


biweight <- function(x,...){
  Nprobes <- dim(x)[1]
  Nchips <- dim(x)[2]
  probes <- as.factor(rep(1:Nprobes,Nchips))
  samps <- as.factor(rep(1:Nchips,rep(Nprobes,Nchips)))
  
  z  <- rlm(as.vector(log2(x))~samps+probes,psi=psi.bisquare,maxit=50,...)
##add se's later
  c(z$coef[1],z$coef[1]+z$coef[2:Nchips],rep(NA,Nchips))
}

avglogpm <- function(x){
  Nprobes <- dim(x)[1]
  Nchips <- dim(x)[2]
  
  probes <- as.factor(rep(1:Nprobes,Nchips))
  samps <- as.factor(rep(1:Nchips,rep(Nprobes,Nchips)))
  
  z  <- lm(as.vector(log2(x))~samps+probes)
  ##add se's later
  c(z$coef[1],z$coef[1]+z$coef[2:Nchips],rep(NA,Nchips))
}




