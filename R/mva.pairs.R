library(modreg)
mva.pairs <- function(x,labels=colnames(x),log.it=TRUE,span=2/3,family.loess="symmetric",digits=3,line.col=2,main="MVA plot",...){
  if(log.it) x <-log2(x)
  J <- dim(x)[2]
  frame()
  old.par <- par(no.readonly = TRUE)
   on.exit(par(old.par))
par(mfrow=c(J,J),mgp=c(0,.2,0),mar=c(1,1,1,1),oma=c(1,1.4,2,1))
  for(j in 1:(J-1)){
    par(mfg=c(j,j));plot(1,1,type="n",xaxt="n",yaxt="n",xlab="",ylab="");text(1,1,labels[j],cex=2)
    for(k in (j+1):J){
      par(mfg=c(j,k))
      yy <- x[,j]-x[,k]
      xx <-(x[,j]+x[,k])/2
      xx <- xx
      yy <- yy
      aux <- loess(yy~xx,degree=1,span=span,family=family.loess)$fitted
      plot(xx,yy,pch=".",xlab="",ylab="",tck=0,...)
      o <- order(xx)
      lines(approx(xx[o],aux[o]),col=line.col)
      par(mfg=c(k,j))
      sigma <- quantile(yy,.75)-quantile(yy,.25)
      txt <- format(c(sigma,0.123456789),digits=digits)
      plot(c(0,1),c(0,1),type="n",ylab="",xlab="",xaxt="n",yaxt="n")
      text(0.5,0.5,txt,cex=2)
    }
  }
  par(mfg=c(J,J));plot(1,1,type="n",xaxt="n",yaxt="n",xlab="",ylab="");
  text(1,1,labels[J],cex=2)
  mtext("A",1,outer=TRUE,cex=1.5)
  mtext("M",2,outer=TRUE,cex=1.5,las=1)
  mtext(main,3,outer=TRUE,cex=1.5)
  invisible()
}
