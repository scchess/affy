.initAffyBatch <- function(where){

  if (debug.affy123) cat("-->initAffyBatch\n")
    
  setClass("AffyBatch",
           representation(intensity="array",
                          cdfName="character",
                          chipNames="character",
                          phenoData="phenoData",
                          nrow="numeric", # dim(intensity) <- c(nrow, ncol, nexp) to go 3-ways
                          ncol="numeric",
                          nexp="numeric",
                          annotation="character",
                          description="characterORmiame",
                          notes="character",
                          history="list"),
           where=where)

  if (debug.affy123) cat("--->accessors\n")
  ## accessors
  setMethod("intensity", signature(object="AffyBatch"),
            function(object) object@intensity,
            where=where)
  
  setReplaceMethod("intensity", signature(object="AffyBatch"),
                   function(object, value){
                     object@intensity <- value
                     object
                   }, where=where)

  setMethod("history", signature(max.show="AffyBatch", reverse="missing"),
            function(max.show) max.show@history,
            where=where)
  
  setReplaceMethod("history", signature(object="AffyBatch"),
                   function(object, value){
                     object@history <- value
                     object
                   }, where=where)

  setMethod("length",signature(x="AffyBatch"),
            function(x) x@nexp,
            where=where)
  
  if (debug.affy123) cat("--->getCdfInfo\n")
  if( !isGeneric("getCdfInfo") )
    setGeneric("getCdfInfo", function(object, ...)
               standardGeneric("getCdfInfo"), where=where)

  setMethod("getCdfInfo", signature("AffyBatch"),
            function (object, what=getOption("BioC")$affy$probesloc.what,
                      where=getOption("BioC")$affy$probesloc.where) {
              ## "what" can be "package", "file" or "environment"
              ## "where" is where it can be found
              
              if (what == "package") {
                loc <- .find.package(object@cdfName, lib.loc=where, quiet=TRUE)
                
                if (identical(loc, character(0)))
                  stop(paste("AffyBatch: Looked for probes information in the package ", object@cdfName, "but could not find it.\n"))
                ##may be an option to try to autoload the package from
                ##the bioconductor website woud be nice here
                
                if (length(loc) > 1)
                  warning(paste("several packages with a matching name. Using the one at", loc[1]))

                do.call("library", list(object@cdfName, lib.loc=dirname(loc[1])))
                return(get(name, envir=as.environment(paste("package:", object@cdfName, sep=""))))
                ##object@cdfInfo <<- get(name, envir=as.environment(paste("package:",name, sep="")))
              }
              
              if (what == "file") {
                cdf <- read.cdffile(file.path(path.expand(where), object@cdfName))
                ## ---> extra paranoia <---
                if (cdf@cdfName != object@cdfName)
                  warning(paste("The CDFALSE file identifies as", cdf@cdfName,
                                "while you probably want", object@cdfName))
                ## ---> end <---
                return(getLocations.Cdf(cdf))
                ##object@cdfInfo <<- getLocations.Cdf(cdf)
                rm(cdf)
                gc() # since cdf can be rather large
              }

              if (what == "environment") {
                return(as.environment(get(object@cdfName, where)))
                ##object@cdfInfo <<- as.environment(get(name, where))
              }
              
            },
            where=where)

  ##geneNames method
  if (debug.affy123) cat("--->geneNames\n")
  if( !isGeneric("geneNames") )
    setGeneric("geneNames", function(object)
               standardGeneric("geneNames"), where=where)
  setMethod("geneNames",signature("AffyBatch"),
            function(object){
              cdf.envir <- getCdfInfo(object)
              return(ls(env=cdf.envir))
            },where=where)

  
  ##chipNames method
  if( !isGeneric("chipNames") )
    setGeneric("chipNames", function(object)
               standardGeneric("chipNames"), where=where)
  setMethod("chipNames","AffyBatch",function(object){
    ##colnames(object@intensity)
    object@chipNames
  },where=where)

  if( !isGeneric("chipNames<-") )
    setGeneric("chipNames<-", function(object)
               standardGeneric("chipNames<-"), where=where)
  setReplaceMethod("chipNames", signature(object="AffyBatch"),
                   function(object, value){
                     ##colnames(intensity(object)) <- value
                     object@chipNames <- value
                     object
                   }, where=where)

  
  ##show method
  if (debug.affy123) cat("--->show\n")
  setMethod("show", "AffyBatch",
            function(object) {
              tmp <- geneNames(object)
              cat("AffyBatch object\n")
              cat("size of arrays=",object@ncol,"x",object@nrow,
                  " features\n",sep="")
              #cat("cdf=",object@cdffile,"\n",sep="")
              cdf.env <- getCdfInfo(object)
              cat("cdf=", object@cdfName,
                  " (", length(ls(env=cdf.env)), " affyids)\n",
                  sep="")
              cat("number of experiments=",length(object),"\n",sep="")
              cat("number of genes=",length(tmp),"\n",sep="")
              cat("annotation=",object@annotation,"\n",sep="")
              cat("notes=",object@notes,"\n",sep="")
              ##cat("sample names=",sampleNames(object),"\n",sep="\n")
            },
            where=where)


  if( !isGeneric("indexProbes") )
    setGeneric("indexProbes", function(object, which, ...)
               standardGeneric("indexProbes"), where=where)

  setMethod("indexProbes", signature("AffyBatch", which="character"),
            function(object, which=c("pm", "mm"), genenames=NULL, xy=FALSE) {
              
              which <- match.arg(which)
              
              i.probes <- match(which, c("pm", "mm")) # i.probes will know if "[,1]" or "[,2]"
              
              envir <- getCdfInfo(object)
              
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
  if (debug.affy123) cat("--->pmindex\n")
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
  if (debug.affy123) cat("--->probeNames\n")
  if( !isGeneric("probeNames") )
    setGeneric("probeNames", function(object, ...)
               standardGeneric("probeNames"), where=where)
  
  setMethod("probeNames","AffyBatch",
  function(object,genenames=NULL,mm=FALSE){
    if(mm) Index <- mmindex(object,genenames)
    else Index <- pmindex(object,genenames)
    reps <- unlist(lapply(Index,length))
    rep(names(Index),reps)
  },where=where)


  if (debug.affy123) cat("--->probes\n")
  if( !isGeneric("probes") )
    setGeneric("probes", function(object, ...)
               standardGeneric("probes"), where=where)
  
  setMethod("probes", signature("AffyBatch"),
            function(object, which=c("pm", "mm"), genenames=NULL, LISTRUE=FALSE){

              which <- match.arg(which)
              
              index <- indexProbes(object, which, genenames)
              
              if(LISTRUE)
                ans <- lapply(index, function(i) object@intensity[i, ])
              else{
                index <- unlist(index)
                ans <- object@intensity[index, ]
                colnames(ans) <- chipNames(object)
                rownames(ans) <- names(index)
              }
              
              return(ans)
            },
            where=where)
  
  ##pm method
  if (debug.affy123) cat("--->pm\n")
  if( !isGeneric("pm") )
    setGeneric("pm", function(object, ...)
               standardGeneric("pm"), where=where)
  
  setMethod("pm","AffyBatch",
            function(object, genenames=NULL, LISTRUE=FALSE) probes(object, "pm", genenames, LISTRUE=LISTRUE),
            where=where
            )
  

  ##mm method
  if( !isGeneric("mm") )
    setGeneric("mm", function(object, ...)
               standardGeneric("mm"), where=where)
  
  setMethod("mm",signature("AffyBatch"),
            function(object, genenames=NULL, LISTRUE=FALSE) probes(object, "mm", genenames, LISTRUE=LISTRUE),
            where=where
            )

  
  if (debug.affy123) cat("--->ppset\n")
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

  if (debug.affy123) cat("--->[[\n")
  setMethod("[[", "AffyBatch",
            function(x, i, j) {
              ##DEBUG: NA ?! watch this in next versions of R
              ## spotsd stuff to be really removed ?
              ##if (is.na.spotsd(x)) {
              mysd <- matrix()
              ##} else {
              ##  mysd <- spotsd(x)[, , i]
              ##}
              oldim <- dim(intensity(x))
              dim(intensity(x)) <- c(x@nrow, x@ncol, x@nexp)
              ##new("Cel", intensity=intensity(x)[, , i], sd=mysd, name=chipNames(x)[i], cdfName=x@cdfName, outliers=outliers(x)[[i]], masks=masks(x)[[i]], history=history(x)[[i]]) ## commented out becuase no 'outliers' or 'masks'.
              new("Cel", intensity=intensity(x)[, , i], sd=mysd, name=chipNames(x)[i], cdfName=x@cdfName, outliers=matrix(), masks=matrix(), history=history(x)[[i]])
              dim(intensity(x)) <- oldim
            },
            where=where)

  setReplaceMethod("[[", "AffyBatch",
                   function(x, i, j, ..., value) {
                     oldim <- dim(intensity(x))
                     dim(intensity(x)) <- c(x@nrow, x@ncol, x@nexp)
                     intensity(x)[, , i] <- intensity(value)
                     ## spotsd ?
                     ##if ((! is.na.spotsd(x)) & (spotsd(value) != c(NA)))
                     ##  spotsd(x)[, , i] <- spotsd(value)
                     x@name[i] <- chipNames(value)
                     
                     if (x@cdfName != value@cdfName)
                       warning("cdfName mismatch !")
                     
                     ##outliers(x)[[i]] <- outliers(value)
                     ##masks(x)[[i]] <- masks(value)
                     history(x)[[i]] <- history(value)
                     dim(intensity(x)) <- oldim
                     return(x)
                   },
                   where=where)
  
  ## normalize.methods
  if (debug.affy123) cat("--->normzalize.methods\n")
  if( !isGeneric("normalize.methods") )
    setGeneric("normalize.methods", function(object)
               standardGeneric("normalize.methods"),
               where=where)
  
  setMethod("normalize.methods", signature(object="AffyBatch"),
            function(object) {
              normalize.AffyBatch.methods
            },
            where=where)

  
  if (! isGeneric("normalize"))
    setGeneric("normalize", function(object, ...) standardGeneric("normalize"),
               where=where)
  
  setMethod("normalize", signature(object="AffyBatch"),
            function(object, f.cdf=NULL, method="quantiles", ...) {
              method <- match.arg(method, normalize.AffyBatch.methods)
              if (is.na(method))
                stop("unknown method")
              method <- paste("normalize.AffyBatch", method, sep=".")
              ## change dimension to re-use Celc.container related code
              oldim <- dim(intensity(object))
              dim(intensity(object)) <- c(object@nrow, object@ncol, object@nexp)
              do.call(method, alist(object, ...))
              dim(intensity(object)) <- oldim
            },
            where=where)

  
  ## expression value computation
  if (debug.affy123) cat("--->normalize\n")
  if( !isGeneric("generatExprSet") )
    setGeneric("generateExprSet", function(x, cdf, method, bg.correct, ...)
               standardGeneric("generateExprSet"),
               where=where)
  
  setMethod("generateExprSet", signature(x="AffyBatch", cdf="missing", method="character", bg.correct="character"),
            function(x, method, bg.correct, ids=NULL, verbose=TRUE, param.bg.correct=list(), param.summary=list()) {
              
              
              ##DEBUG: check the existence of method and bg.correct HERE !
              n <- length(x)

              ## if NULL compute for all
              if (is.null(ids))
                ids <- geneNames(x)
              
              m <- length(ids)
              
              idsi <- match(ids, geneNames(x))
              
              
              ## cheap trick to (try to) save time
              c.pps <- new("PPSet.container",
                           pmProbes=matrix(),
                           mmProbes=matrix())

              
              ## matrix to hold expression values
              exp.mat <- matrix(NA, m, n)

              
              ##if (verbose) cat("getting the locations for all affyIDs.....")
              ##all.l.pm <- .Call("getallLocations", as.integer(cdf@name), as.integer(dim(cdf@name)),
              ##                  as.integer(atom(cdf)), as.integer(pmormm(cdf)),
              ##                  as.integer(max(cdf@name, na.rm=TRUE)+1) )
              ##all.l.mm <- .Call("getallLocations", as.integer(cdf@name), as.integer(dim(cdf@name)),
              ##                  as.integer(atom(cdf)), as.integer(! pmormm(cdf)),
              ##                  as.integer(max(cdf@name, na.rm=TRUE)+1) )

              if (verobse) cat(".....done.\n")
              
              ##DEBUG: hackish (put global adjsutment names below
              if (bg.correct %in% c("bg.correct.rma")) {
                if (verbose) cat("computing parameters for global background adjustement.....")
                all.l.pm.mat <- unlist(lapply(multiget(ids, env=getCdfInfo(x)),  function(x) if (ncol(x) == 2) x[,1]))
                all.l.mm.mat <- unlist(lapply(multiget(ids, env=getCdfInfo(x)),  function(x) if (ncol(x) == 2) x[,2]))
                ##all.l.pm.mat <- cbind(unlist(lapply(all.l.pm, function(x) if (ncol(x) == 2) x[,1])),
                ##                      unlist(lapply(all.l.pm, function(x) if (ncol(x) == 2) x[,2])))
                ##all.l.mm.mat <- cbind(unlist(lapply(all.l.mm, function(x) if (ncol(x) == 2) x[,1])),
                ##                      unlist(lapply(all.l.mm, function(x) if (ncol(x) == 2) x[,2])))
                all.param <- lapply(seq(1:n), function(i) {
                  notNA <- !(is.na(intensity(x)[, i][all.l.pm.mat]) | is.na(intensity(x)[, i][all.l.mm.mat]))
                  bg.parameters(intensity(x)[, i][notNA], intensity(x)[, i][notNA])
                })
                param.bg.correct$all.param <- all.param
                if (verbose) cat(".....done.\n")
              }
              
             if (verbose) {
                cat(m,"ids to be processed\n")
                countprogress <- 0
              }
              
              ## loop over the ids
              mycall <- as.call(c(getMethod("express.summary.stat", signature=c("PPSet.container","character","character")),
                                  list(c.pps, method=method, bg.correct=bg.correct, param.bg.correct=param.bg.correct, param.method=param.summary)))

              options(show.error.messages = FALSE)
              on.exit(options(show.error.messages = TRUE))
              
              for (i in seq(along=ids)) {
                
                id <- ids[i]

                ##cat(i,"--")
                if (verbose) {
                  #cat(id, "\n")
                  if ( round(m/10) == countprogress) {
                    cat(".")
                    countprogress <- 0
                  }
                  else
                    countprogress <- countprogress + 1
                }
                ## locations for an id
                ##l.pm <- locate.name(ids[id], cdf, type="pm")
                ##l.mm <- locate.name(ids[id], cdf, type="mm")
                loc <- get(id ,envir=getCdfInfo(x))
                l.pm <- loc[, 1]
                if (ncol == 2)
                  l.mm <- loc[ ,2]
                else
                  l.mm <- NA
                
                ## fill the PPSet.container
                ##c.pps@pmProbes <- matrix(NA, nrow=length(l.pm), ncol=n)
                ##c.pps@mmProbes <- matrix(NA, nrow=length(l.pm), ncol=n)

                np <- length(l.pm)
                
                ##names are skipped

                warning("indexing of probes not yet checked (possibly wrong)")
                c.pps@pmProbes <- matrix(intensity(x)[l.pm, ],
                                         np, n, byrow=TRUE)
                c.pps@mmProbes <- matrix(intensity(x)[l.mm, ],
                                         np, n, byrow=TRUE)
                
                ##c.pps@pmProbes <- matrix(intensity(x)[cbind(matrix(rep(l.pm, rep(n, np*2)), nrow=np*n, ncol=2, byrow=FALSE),
                ##                                            rep(1:n, np))],
                ##                         np, n, byrow=TRUE)
                ##c.pps@mmProbes <- matrix(intensity(x)[cbind(matrix(rep(l.mm, rep(n, np*2)), nrow=np*n, ncol=2, byrow=FALSE),
                ##                                            rep(1:n, np))],
                ##                         np, n, byrow=TRUE)
                
                ## generate expression values
                ## (wrapped in a sort of try/catch)
                mycall[[2]] <- c.pps
                ev <- try(eval(mycall))
                
                if (! inherits(ev,"try-error")) {
                  exp.mat[i, ] <- ev
                } else {
                  warning(paste("Error with affyid:", name.levels(cdf)[id]))
                }
                ## no need for an 'else' branching since exp.mat was initialized with NA
                
              }
              
              options(show.error.messages = TRUE)
              on.exit(NULL)
              
              if (verbose) cat("\n")

              ## instance exprSet
              ##if (verbose) cat("instancianting an exprSet.....")
              dimnames(exp.mat) <- list(ids, deparse(x))
              eset <- new("exprSet", exprs=exp.mat, se.exprs=matrix())
              ##if (verbose) cat(".....done.\n")
              
              return(eset)
            },
            where=where)
  
  ## use [[ and image instead !
  
 #  if( !isGeneric("image") )
#     setGeneric("image",where=where)
#   setMethod("image",signature(x="AffyBatch"),
#             function(x, transfo=log, ...){
#               scn <- prod(par("mfrow"))
#               ask <- dev.interactive()
#               which.plot <- 0
#               for(i in 1:length(sampleNames(x))){
#                 which.plot <- which.plot+1;
#                 if(trunc((which.plot-1)/scn)==(which.plot-1)/scn && which.plot>1 && ask)  par(ask=TRUE)
#                 image(AffyBatch[[i]], ...)
#                 par(ask=FALSE)}
#             },where=where)
  

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
#    if(LISTRUE)
#      ans <- lapply(Index,function(i) object@intensity[i,])
#    else{
#      Index2 <- unlist(Index)
#      ans <- object@intensity[Index2,]
#      colnames(ans) <- sampleNames(object)
#      rownames(ans) <- names(Index2)
#    }
#    return(ans)
#  },where=where)

#  setMethod("mm","AffyBatch",function(object,genenames=NULL,LISTRUE=FALSE){
#    Index <- mmindex(object,genenames)
#    if(LISTRUE)
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





