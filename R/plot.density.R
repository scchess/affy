plot.density <- function(x,col=rainbow(dim(x)[2]),ylab="density",xlab="x",lwd=rep(2,dim(x)[2]),lty=rep(1,dim(x)[2]),xlim=range(x),ylim=c(0,1),...){

  max.y <- function(dens){
	max(dens$y)
  }
	
  x.density <- apply(x,2,density)
  
  if ((min(ylim) == 0) & (max(ylim) == 1)){
  	max.ylim <-lapply(x.density,max.y)
	ylim <- c(0,max(unlist(max.ylim)))
  }

 plot(c(0,0),c(0,0),type="n",ylab=ylab,xlab=xlab,xlim=xlim,ylim=ylim,...)
 for (i in 1:length(x.density)){
  lines(x.density[[i]],col=col[i],lwd=lwd[i],lty=lty[i],...)
 }
}
 

plot.density.AffyBatch <- function(x,col=rainbow(length(x)),log=TRUE,
                                   which=c("pm","mm","both"),
                                   ylab="density",xlab=NULL,
                                   lwd=rep(2,length(x)),
                                   lty=rep(1,length(x)),
                                   xlim=NULL,
                                   ylim=c(0,1),...){
  Index <- unlist(indexProbes(x,which=which))
  x <- intensity(x)[Index,]
  if(log){
    x <- log2(x)
    if(is.null(xlab)) xlab <- "log intensity"
  }
  else  if(is.null(xlab)) xlab <- "intensity"

  if(is.null(xlim)) xlim <- range(x)
  
  plot.density(x,col=col,ylab=ylab,xlab=xlab,
               lwd=lwd,
               lty=lty,
               xlim=xlim,
               ylim=ylim,...)
}
