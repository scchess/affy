####These functions take pm and mm vectors and return background corrected
####pms and mms. they will be then converted to expression by summary
####functions
bg.parameters <- function(pm, mm, n.pts=2^14){

  max.density <- function(x,n.pts){
    aux <- density(x,kernel="epanechnikov",n=n.pts)
    aux$x[order(-aux$y)[1]] 
  }
  
  mmbg <- max.density(mm,n.pts)
  pmbg <- max.density(pm,n.pts)
 
  bg.data <- mm[mm < mmbg]
  bg.data <- bg.data-mmbg
  bgsd <- sqrt(sum(bg.data^2)/(length(bg.data)-1))*sqrt(2)/.85
 
  sig.data <- pm[pm > pmbg]
  sig.data <- sig.data-pmbg
 
  alpha <- 1/mean(sig.data)
  mubg <- mmbg
  list(alpha=alpha,mu=mubg,sigma=bgsd)  
}

bg.adjust <- function(pm, mm, n.pts=2^14){
  param <- bg.parameters(pm, mm, n.pts)
  b <- param$sigma
  a <- pm - param$mu - param$alpha*b^2
  a + b*((1./sqrt(2*pi))*exp((-1./2.)*((a/b)^2)))/pnorm(a/b)
}

subtractmm <- function(pm, mm){
  return(pm - mm)
}

bgc <- function(object,bg=bg.adjust){
  n <- ncol(pm(object))
  devnull <- lapply(1:n,
                    function(i) {
                      pm(object)[,i] <<- bg(pm(object)[,i], mm(object)[,i]);
                      NULL}
                    )
  object
}

##DEBUG: merge with what is found before...
bg.correct.subtractmm <- function(pm, mm){
  return(pm-mm)
}

bg.correct.pmonly <- function(pm, mm) {
  return(pm)
}

bg.correct.adjust <- function(pm ,mm, all.param){
  r <- matrix(NA, nrow(pm), ncol(pm))
  #cat(str(r))
  for (i in 1:ncol(r)) {
    b <- all.param[[i]]$sigma
    a <- pm[,i] - all.param[[i]]$mu - all.param[[i]]$alpha*b^2
    #cat(str(a + b*((1./sqrt(2*pi))*exp((-1./2.)*((a/b)^2)))/pnorm(a/b)))
    r[, i] <- a + b*((1./sqrt(2*pi))*exp((-1./2.)*((a/b)^2)))/pnorm(a/b)
  }
  return(r)
}
