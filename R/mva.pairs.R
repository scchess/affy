###
###
### Code for M and MvA plots
###
### Mar 6, 2004 - added the generic Mbox. It performs
###               the equivalent of Mbox in affyPLM
###               added a generic MAplot. Similar
###               functionality is implemented in
###               affyPLM
###               a function ma.plot now does the actual plotting
###               for mva.pairs
###
### Aug 23, 2004 - change the placement location of statistics in
###                ma.plot

ma.plot <- function(A,M,subset=sample(1:length(M),min(c(10000, length(M)))),show.statistics=TRUE,span=2/3,family.loess="gaussian",cex=2,...){

  fn.call <- list(...)

  sigma <- IQR(M)
  mean <- median(M)
  if (!is.element("ylim",names(fn.call))){
    yloc <- max(M)
  } else {
    yloc <- max(fn.call$ylim)
  }
  if (!is.element("xlim",names(fn.call))){
    xloc <- max(A)
  } else {
    yloc <- max(fn.call$xlim)
  }
  
  aux <- loess(M[subset]~A[subset],degree=1,span=span,family=family.loess)$fitted
  
  plot(A,M,...)
  o <- order(A[subset])
  A <- A[subset][o]
  M <- aux[o]
  o <-which(!duplicated(A))
  lines(approx(A[o],M[o]),col="red")
  abline(0,0,col="blue")

  # write IQR and Median on to plot
  if (show.statistics){
    txt <- format(sigma,digits=3)
    txt2 <- format(mean,digits=3)
    text(xloc ,yloc,paste(paste("Median:",txt2),paste("IQR:",txt),sep="\n"),cex=cex,adj=c(1,1))
  }

  
}



mva.pairs <- function(x,labels=colnames(x),log.it=TRUE,span=2/3,family.loess="gaussian",digits=3,line.col=2,main="MVA plot",...){
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
      sigma <- IQR(yy)
      mean <- median(yy)
      subset<-sample(1:length(x),min(c(10000, length(x))))
      ma.plot(xx,yy,tck=0,subset=subset,show.statistics=FALSE,pch=".",xlab="",ylab="",tck=0,...)
      par(mfg=c(k,j))
      #sigma <- IQR(yy)
      txt <- format(sigma,digits=digits)
      txt2 <- format(mean,digits=digits)
      plot(c(0,1),c(0,1),type="n",ylab="",xlab="",xaxt="n",yaxt="n")
      text(0.5,0.5,paste(paste("Median:",txt2),paste("IQR:",txt),sep="\n"),cex=2)
    }
  }
  par(mfg=c(J,J));plot(1,1,type="n",xaxt="n",yaxt="n",xlab="",ylab="");
  text(1,1,labels[J],cex=2)
  mtext("A",1,outer=TRUE,cex=1.5)
  mtext("M",2,outer=TRUE,cex=1.5,las=1)
  mtext(main,3,outer=TRUE,cex=1.5)
  invisible()
}









if (!isGeneric("Mbox"))
  setGeneric("Mbox",function(object,...)
             standardGeneric("Mbox"))


setMethod("Mbox",signature("AffyBatch"),
          function(object,log=TRUE,type=c("both","pm","mm"),...){
             type <- match.arg(type)
             if (type == "both"){
              pms <- unlist(indexProbes(object, "both"))
            } else if (type == "pm"){
              pms <- unlist(pmindex(object))
            } else if (type == "mm"){
              mms <- unlist(mmindex(object))
            }
            if(log){
              x <- log2(intensity(object)[pms, ])
            } else {
              x <- intensity(object)[pms, ]
            }
            medianchip <- apply(x, 1, median)
            M <- sweep(x,1,medianchip,FUN='-')
            boxplot(data.frame(M),...)
          })

if (!isGeneric("MAplot"))
  setGeneric("MAplot",function(object,...)
             standardGeneric("MAplot"))


setMethod("MAplot",signature("AffyBatch"),
          function(object,log=TRUE,type=c("both","pm","mm"),ref=NULL,subset=NULL,which=NULL,...){
            type <- match.arg(type)
            if (type == "both"){
              pms <- unlist(indexProbes(object, "both"))
            } else if (type == "pm"){
              pms <- unlist(pmindex(object))
            } else if (type == "mm"){
              pms <- unlist(mmindex(object))
            }
            if(log){
              x <- log2(intensity(object)[pms, ])
            } else {
              x <- intensity(object)[pms, ]
            }

            if (is.null(which)){
              which <- 1:dim(object@exprs)[2]
            }
            

            if (is.null(subset)){
              if (is.null(ref)){
                medianchip <- apply(x, 1, median)
              } else {
                medianchip <- x[,ref]
              }
              M <- sweep(x,1,medianchip,FUN='-')
              A <- 1/2*sweep(x,1,medianchip,FUN='+')
              if (is.null(ref)){
                for (i in which){
                  title <- paste(sampleNames(object)[i],"vs pseudo-median reference chip")
                  ma.plot(A[,i],M[,i],main=title,xlab="A",ylab="M",pch='.',...)
                }
              } else {
                for (i in which){
                  if (which != ref){
                    title <- paste(sampleNames(object)[i],"vs",sampleNames(object)[ref])
                    ma.plot(A[,i],M[,i],main=title,xlab="A",ylab="M",pch='.',...)
                  }
                }
              }
            } else {
              if (is.null(ref)){
                medianchip <- apply(x[,subset], 1, median)
              } else {
                if (is.element(ref,subset)){
                  medianchip <- x[,ref]
                } else {
                  stop("Ref ",ref, "is not part of the subset")
                }
              }
              if (!all(is.element(which,subset))){
                stop("Specified arrays not part of subset")
              }
              M <- sweep(x,1,medianchip,FUN='-')
              A <- 1/2*sweep(x,1,medianchip,FUN='+')
              if (is.null(ref)){
                for (i in which){
                  title <- paste(sampleNames(object)[i],"vs pseudo-median reference chip")
                  ma.plot(A[,i],M[,i],main=title,xlab="A",ylab="M",pch='.',...)
                }
              } else {
                for (i in which){
                  if (i != ref){
                    title <- paste(sampleNames(object)[i],"vs",sampleNames(object)[ref])
                    ma.plot(A[,i],M[,i],main=title,xlab="A",ylab="M",pch='.',...)
                  }
                }
              }
            }
          })
