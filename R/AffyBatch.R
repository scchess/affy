.initAffyBatch <- function(where){

  setClass("AffyBatch",
           representation(intensity="array",
                          cdfName="character",
                          chipNames="character",
                          phenoData="phenoData",
                          nrow="numeric", # dim(intensity) == c(nrow, ncol, nexp)
                          ncol="numeric",
                          nexp="numeric",
                          annotation="character",
                          description="characterORmiame",
                          notes="character",
                          history="list"),
           where=where)
  
  setMethod("intensity", signature(object="AffyBatch"),
            function(object) object@intensity,
            where=where)
  setReplaceMethod("intensity", signature(object="AffyBatch"),
                   function(object, value){
                     object@intensity <- value
                     object
                   }, where=where)

  if( !isGeneric("getCdfInfo") )
    setGeneric("getCdfInfo", function(object, name, ...)
               standardGeneric("getCdfInfo"), where=where)

  setMethod("getCdfInfo", signature("AffyBatch", "character"),
            function (object, name, what="package", where=NULL, compress=FALSE) {
              ## "what" can be "package", "file" or "environment"
              ## "where" is where it can be found
              ## "compress" only makes sense when reading from a CDF file
              
              if (what == "package") {
                loc <- .find.package(name, lib.loc=where, quiet=TRUE)
                
                if (identical(loc, character(0)))
                  stop(paste("No package like", name, "could be found.\n"))
                ##may be an option to try to autoload the package from
                ##the bioconductor website woud be nice here
                
                if (length(loc) > 1)
                  warning(paste("several packages with a matching name. Using the one at", loc[1]))

                do.call("library", list(name, lib.loc=dirname(loc[1])))
                return(get(name, envir=as.environment(paste("package:",name, sep=""))))
                ##object@cdfInfo <<- get(name, envir=as.environment(paste("package:",name, sep="")))
                ##invisible(TRUE)
              }
              
              if (what == "file") {
                cdf <- read.cdffile(file.path(path.expand(where), name), compress=compress)
                ## ---> extra paranoia <---
                if (cdf@cdfName != object@cdfName)
                  warning(paste("The CDF file identifies as", cdf@cdfName,
                                "while you probably want", object@cdfName))
                ## ---> end <---
                return(getLocations.Cdf(cdf))
                ##object@cdfInfo <<- getLocations.Cdf(cdf)
                rm(cdf)
                gc() # since cdf can be rather large
                ##invisible(TRUE)
              }

              if (what == "environment") {
                return(as.environment(get(name, where)))
                ##object@cdfInfo <<- as.environment(get(name, where))
                ##invisible(TRUE)
              }
              
            },
            where=where)

#  ##loadcdf method
#  if( !isGeneric("loadcdf") )
#    setGeneric("loadcdf", function(object,...)
#               standardGeneric("loadcdf"), where=where)
                                        #  ##loadcdf method
#  setMethod("loadcdf","AffyBatch",function(object,lib.loc=NULL){
#    ##try to load library.
#    ##we know which one becuase its in object@cdffile
#    options(show.error.messages = FALSE)
#    tmp <- try(do.call("library", list(package=object@cdffile,lib.loc=lib.loc)))
#    options(show.error.messages = TRUE)
    
#    ##if not there install and download, well for now just give error
#    if(inherits(tmp, "try-error")){
#      stop(paste("You need to dowload and install the",object@cdffile,"library.\n"))
#    }
#    ##return the evnironmet
#    get(object@cdffile,envir=get(object@cdffile))
#  },where=where)


  ##geneNames method
  if( !isGeneric("geneNames") )
    setGeneric("geneNames", function(object)
               standardGeneric("geneNames"), where=where)
  setMethod("geneNames","AffyBatch",function(object){
    return(ls(env=object@cdfInfo))
  },where=where)


  ##sampleNames method
  if( !isGeneric("sampleNames") )
    setGeneric("sampleNames", function(object)
               standardGeneric("sampleNames"), where=where)
  setMethod("sampleNames","AffyBatch",function(object){
    colnames(object@intensity)
  },where=where)


  
  ##show method
  setMethod("show", "AffyBatch",
            function(object) {
              tmp <- geneNames(object)
              cat("AffyBatch object\n")
              cat("size of arrays=",object@ncol,"x",object@nrow,
                  " features\n",sep="")
              #cat("cdf=",object@cdffile,"\n",sep="")
              cat("cdf=", object@cdfName,
                  " (", length(ls(env=object@cdfInfo)), " affyids)\n",
                  sep="")
              cat("number of samples=",length(sampleNames(object)),"\n",sep="")
              cat("number of genes=",length(tmp),"\n",sep="")
              cat("annotation=",object@annotation,"\n",sep="")
              cat("notes=",object@notes,"\n",sep="")
              cat("sample names=",sampleNames(object),"\n",sep="\n")
            },
            where=where)


  if( !isGeneric("indexProbes") )
    setGeneric("indexProbes", function(object, which, ...)
               standardGeneric("indexProbes"), where=where)

  setMethod("indexProbes", signature("AffyBatch", which="character"),
            function(object, which=c("pm", "mm"), genenames=NULL, xy=FALSE) {
              
              which <- match.arg(which)
              
              i.probes <- match(which, c("pm", "mm")) # i.probes will know if "[,1]" or "[,2]"
              
              envir <- object@cdfInfo
              
              if(is.null(genenames)) 
                genenames <- ls(envir )
              
              ## shorter code, using the features of multiget
              ## (eventually more readable too)
              ## note: genenames could be confusing (the same gene can be
              ## found in several affyid (ex: the 3' and 5' controls)
              
              ans <-  multiget(genenames, pos, envir, iffail=NA)

              ## this kind of thing could be included in 'multiget' as
              ## and extra feature. A function could be specified to
              ## process what is 'multiget' on the fly
              for (i in seq(along=ans)) {
                
                if ( is.na(ans[[i]]) )
                  next
                
                tmp <- ans[[i]][, i.probes]
                
                if (xy) {
                  x <- tmp %% nrow(object)
                  y <- tmp %/% nrow(object) + 1
                  ans[[i]] <- cbind(x, y)
                }
                
                ans[[i]] <- tmp
              }
              
              return(ans)
            },
            where=where)
  
  
  ##pmindex method
  if( !isGeneric("pmindex") )
    setGeneric("pmindex", function(object,...)
               standardGeneric("pmindex"), where=where)

  ##wrapper
  setMethod("pmindex", "AffyBatch",
            function(object,genenames=NULL,xy=FALSE) 
            indexProbes(object, "pm", genenames=genenames, xy=xy),
            where=where
            )
            
  ##mmindex method
            if( !isGeneric("mmindex") )
            setGeneric("mmindex", function(object,...)
               standardGeneric("mmindex"), where=where)
  
  ##wrapper
  setMethod("mmindex", "AffyBatch",
            function(object,genenames=NULL,xy=FALSE) 
            indexProbes(object, "mm", genenames=genenames, xy=xy),
            where=where
            )                        

  
  ##probeNames method
  if( !isGeneric("probeNames") )
    setGeneric("probeNames", function(object, ...)
               standardGeneric("probeNames"), where=where)
  
  setMethod("probeNames","AffyBatch",
  function(object,genenames=NULL,mm=F){
    if(mm) Index <- mmindex(object,genenames)
    else Index <- pmindex(object,genenames)
    reps <- unlist(lapply(Index,length))
    rep(names(Index),reps)
  },where=where)


  if( !isGeneric("probes") )
    setGeneric("probes", function(object, ...)
               standardGeneric("probes"), where=where)

  setMethod("probes", signature("AffyBatch", which="character"),
            function(object, which=c("pm", "mm"), genenames=NULL, LIST=F){

              which <- match.arg(which)
              
              index <- indexProbes(object, which, genenames)
              
              if(LIST)
                ans <- lapply(index, function(i) object@intensity[i, ])
              else{
                index <- unlist(index)
                ans <- object@intensity[index, ]
                colnames(ans) <- sampleNames(object)
                rownames(ans) <- names(index)
              }
              
              return(ans)
            },
            where=where)
  
  cat("--->")
  ##pm method
  if( !isGeneric("pm") )
    setGeneric("pm", function(object, ...)
               standardGeneric("pm"), where=where)
  
  setMethod("pm","AffyBatch",
            function(object, genenames=NULL, LIST=F) probes(object, "pm", genenames, LIST=LIST),
            where=where
            )
  

  ##mm method
  if( !isGeneric("mm") )
    setGeneric("mm", function(object, ...)
               standardGeneric("mm"), where=where)
  
  setMethod("mm",signature("AffyBatch"),
            function(object, genenames=NULL, LIST=F) probes(object, "mm", genenames, LIST=LIST),
            where=where
            )

  
  if( !isGeneric("ppset") )
    setGeneric("ppset", function(object, ...)
               standardGeneric("ppset"), where=where)
  
  setMethod("ppset", "AffyBatch", function(object, genenames=NULL){
    envir <- loadcdf(object)
    if(is.null(genenames))
      genenames <- ls(envir)

    p.pps <- vector("list", length(genenames))
    names(p.pps) <- genenames
    
    for (i in seq(along=genenames)) {
      
      i.pm <- indexProbes(object, "pm", genenames[i])[[1]]
      if (is.na(i.pm))
        intensity.pm <- NA
      else
        intensity.pm <- object(intensity[i.pm, ])
      
      i.mm <- indexProbes(object, "mm", genenames[i])[[1]]
      if (is.na(i.mm))
        intensity.mm <- NA
      else
        intensity.mm <- object(intensity[i.pm, ])
      
      intensity.mm <- object(intensity[i.pm, ])
      p.pps[[i]] <- new("PPSetBatch",
                        pmProbes = intensity.pm,
                        mmProbes = intensity.mm)
    }
    
    return(p.pps)
    
  },where=where)


  setMethod("[[", "AffyBatch",
            function(x, i, j) {
              ##DEBUG: NA ?! watch this in next versions of R
              ## spotsd stuff to be really removed ?
              ##if (is.na.spotsd(x)) {
              mysd <- matrix()
              ##} else {
              ##  mysd <- spotsd(x)[, , i]
              ##}
              new("Cel", intensity=intensity(x)[, , i], sd=mysd, name=x@name[i], cdfName=x@cdfName, outliers=outliers(x)[[i]], masks=masks(x)[[i]], history=history(x)[[i]])
            },
            where=where)

  setReplaceMethod("[[", "AffyBatch",
                   function(x, i, j, ..., value) {
                     intensity(x)[, , i] <- intensity(value)
                     ## spotsd ?
                     ##if ((! is.na.spotsd(x)) & (spotsd(value) != c(NA)))
                     ##  spotsd(x)[, , i] <- spotsd(value)
                     x@name[i] <- chipNames(value)
                     
                     if (x@cdfName != value@cdfName)
                       warning("cdfName mismatch !")
                     
                     outliers(x)[[i]] <- outliers(value)
                     masks(x)[[i]] <- masks(value)
                     history(x)[[i]] <- history(value)
                     return(x)
                   },
                   where=where)
  
  
  if( !isGeneric("image") )
    setGeneric("image",where=where)
  setMethod("image",signature(x="AffyBatch"),
            function(x, transfo=log, ...){
              scn <- prod(par("mfrow"))
              ask <- dev.interactive()
              which.plot <- 0
              for(i in 1:length(sampleNames(x))){
                which.plot <- which.plot+1;
                if(trunc((which.plot-1)/scn)==(which.plot-1)/scn && which.plot>1 && ask)  par(ask=T)
                image(AffyBatch[[i]], ...)
                par(ask=F)}
            },where=where)
  

 #  setMethod("pmindex", "AffyBatch", function(object,genenames=NULL,xy=FALSE){
#    envir <- loadcdf(object)
    
#    if(is.null(genenames)) 
#      genenames <- ls(envir )
     
#    ##this is multiget except we return only pm indeces and the xy condition
#    lengenenames <- length(genenames)
#    ans <- vector("list", length=lengenenames)
    
#    if( ! is.environment(envir) )
#      stop("The cdffile slot does not define an environment")
    
#    opt.err <- options(show.error.messages)
#    options(show.error.messages = FALSE)
#    on.exit(options(show.error.messages = opt.err))
    
#    if (xy)
#      for(i in 1:lengenenames){
#        tmp <-  try(get(genenames[i],pos,envir, "any", TRUE)[,1])
#        y <- floor(tmp/object@nrow)
#        ans[[i]] <- cbind(x=tmp-y*object@nrow-1,y)
#      }
#    else
#      for(i in 1:lengenenames)
#        ans[[i]] <- try(get(genenames[i],pos,envir, "any", TRUE)[,1])
    
#    options(show.error.messages = opt.err)
#    on.exit(NULL)
    
#    failfun <- function(x) {
#      cx <- class(x)
#      if( !is.null(cx) && cx == "try-error")
#        TRUE
#      else
#        FALSE
#    }
#    failed <- sapply(ans, failfun)
#    ans[failed] <- NA
    
#    names(ans) <- genenames
#    ans
#  },where=where) 

#  setMethod("mmindex","AffyBatch",function(object,genenames=NULL,xy=FALSE){
#    envir <- loadcdf(object)
#    if(is.null(genenames))
#      genenames <- ls(envir)
    
#    ##this is multiget except we return only mm indeces and the xy condition
#    lengenenames <- length(genenames)
#    ans <- vector("list", length=lengenenames)
#    if( ! is.environment(envir) )
#      stop("The cdffile slot does not define an environment")
#    options(show.error.messages = FALSE)
#    on.exit(options(show.error.messages = TRUE))
#    if(xy)
#      for(i in 1:lengenenames){
#        tmp <-  try(get(genenames[i],pos,envir, "any", TRUE)[,2])
#        y <- floor(tmp/object@nrow)
#        ans[[i]] <- cbind(x=tmp-y*object@nrow-1,y)
#      }
#    else
#      for(i in 1:lengenenames)
#        ans[[i]] <- try(get(genenames[i],pos,envir, "any", TRUE)[,2])
#    options(show.error.messages = TRUE)
#    on.exit(NULL)
    
#    failfun <- function(x) {
#      cx <- class(x)
#      if( !is.null(cx) && cx == "try-error")
#        TRUE
#      else
#        FALSE
#    }
#    failed <- sapply(ans, failfun)
#    ans[failed] <- NA
    
#    names(ans) <- genenames
#    ans
#  },where=where)
  
#  setMethod("pm","AffyBatch",function(object, genenames=NULL){
#    Index <- pmindex(object,genenames)
#    if(LIST)
#      ans <- lapply(Index,function(i) object@intensity[i,])
#    else{
#      Index2 <- unlist(Index)
#      ans <- object@intensity[Index2,]
#      colnames(ans) <- sampleNames(object)
#      rownames(ans) <- names(Index2)
#    }
#    return(ans)
#  },where=where)

#  setMethod("mm","AffyBatch",function(object,genenames=NULL,LIST=F){
#    Index <- mmindex(object,genenames)
#    if(LIST)
#      ans <- lapply(Index,function(i) object@intensity[i,])
#    else{
#      Index2 <- unlist(Index)
#      ans <- object@intensity[Index2,]
#      colnames(ans) <- sampleNames(object)
#      rownames(ans) <- names(Index2)
#    }
#    return(ans)
#  },where=where)  
  
}





