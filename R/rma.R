##this function does the expression for our method: RMA. not flexible
rma <- function(object,subset=NULL, phenodata=NULL, annotation=NULL,
                description=NULL, notes=NULL, verbose=TRUE, ...){

  
  bg.parameters2 <-  function(pm,n.pts=2^14){
    max.density <- function(x,n.pts){
      aux <- density(x,kernel="epanechnikov",n=n.pts)
      aux$x[order(-aux$y)[1]] 
    }
    pmbg <- max.density(pm,n.pts) ##Log helps detect mode
    
    bg.data <- pm[pm < pmbg]
    ##do it again to really get the mode
    pmbg <- max.density(bg.data,n.pts) 
    bg.data <- pm[pm < pmbg]
    bg.data <- bg.data - pmbg
  
    bgsd <- sqrt(sum(bg.data^2)/(length(bg.data)-1))*sqrt(2)#/.85
    
    sig.data <- pm[pm > pmbg]
    sig.data <- sig.data-pmbg
    
    expmean <- max.density(sig.data,n.pts)
    alpha <- 1/expmean
    mubg <- pmbg
    list(alpha=alpha,mu=mubg,sigma=bgsd)  
  }
  
  bg.adjust2 <- function(pm,n.pts=2^14){
    param <- bg.parameters2(pm,n.pts)
    b <- param$sigma
    a <- pm - param$mu - param$alpha*b^2
    a + b*((1./sqrt(2*pi))*exp((-1./2.)*((a/b)^2)))/pnorm(a/b)
  }
  
  if(verbose) cat("Background correcting\n")
  x <- apply(pm(object),2,bg.adjust2)
  if(verbose) cat("Normalizing Data\n")
  x <- normalize.quantiles(x)
  
  if(verbose) cat("Preparing Data\n")
  if(is.null(subset)){
    data.lst <-split(as.vector(x),factor(rep(probeNames(object),nprobes(object)/dim(x)[1]*dim(x)[2])))
  }
  else{
    if(is.numeric(subset)) subset <- geneNames(object)[subset]
    Index <- which(probeNames(object)%in%subset)
    x <- x[Index,]
    Names <- probeNames(object)[Index]
    data.lst <- split(as.vector(x), factor(rep(Names,length(Names)/dim(x)[1]*dim(x)[2] )))
  }
  data.lst<-lapply(data.lst,matrix,ncol=nchips(object))
  if(verbose) cat("Computing expression. This may take a while.\n")
  e <- t(sapply(data.lst,medianpolish,...))

  if(is.null(phenodata)) phenodata <- phenoData(object)
  if(is.null(annotation)) annotation <- annotation(object)
  if(is.null(description)) description <- description(object)
  if(is.null(notes)) notes <- notes(object)

  exprs <- e[,1:nchips(object)]
  colnames(exprs) <- sampleNames(object)
  se.exprs <- e[,(nchips(object)+1):(2*nchips(object))] 
  colnames(se.exprs) <- sampleNames(object)

  new("exprSet",
      exprs=exprs,
      se.exprs=se.exprs,
      phenoData=phenodata,
      annotation=annotation,
      description=description,
      notes=notes)
}

  


