.initAffyBatch <- function(where){
  
  if (debug.affy123) cat("-->initAffyBatch\n")
  
  setClass("AffyBatch", ##keep it very simple and like exprSet
           representation(cdfName="character",
                          nrow="numeric",
                          ncol="numeric"),
           prototype=list(exprs=matrix(nr=0,nc=0),
             se.exprs = matrix(nr=0,nc=0),
             description=new("MIAME"),
             annotation="",
             notes="",
             cdfName="",
             nrow=0,
             ncol=0), contains="exprSet", where=where)

  if (debug.affy123) cat("--->accessors\n")

#######################################################
### accessors
#######################################################
  
  ##intensity
  setMethod("intensity", signature(object="AffyBatch"),
            function(object) object@exprs,
            where=where)
  
  
  setReplaceMethod("intensity", signature(object="AffyBatch"),
                   function(object, value){
                     object@exprs <- value
                     colnames(object@exprs) <- sampleNames(object)
                     object
                   }, where=where)

  ##for now, there is no accessor for se.exprs. we could use this to store
  ##sd, but if no one uses it... why do it
  
  setMethod("length",signature(x="AffyBatch"),
            function(x) ncol(x@exprs), ##RI: assumes matrices
            where=where)
  
  ##sample Names now comes from Biobase
  if (debug.affy123) cat("--->getCdfInfo\n")
  if( !isGeneric("getCdfInfo") )
    setGeneric("getCdfInfo", function(object, ...)
               standardGeneric("getCdfInfo"), where=where)
  
  setMethod("getCdfInfo", signature("AffyBatch"),
            function (object, how=getOption("BioC")$affy$probesloc) {
              ## "how" is a list. each element of the list must have an element
              ## tagged "what" and an element tagged "where"
              ## "what" can be "data", "package", "file" or "environment"
              ## "where" is where it can be found

              cdfname <- cleancdfname(object@cdfName)
              
              for (i in 1:length(how)) {

                what <- how[[i]]$what
                where <- how[[i]]$where

                if (what == "data") {
                  ##if we can get it from data dir. otherwise load package
                  if(cdfname%in%data(package=affy)$results[,3]){
                    ##RI: package="affy" doesnt work it has to be package=affy
                    ##    fix if you can
                    where <- as.environment(match(paste("package:", where, sep = ""),search()))
                    if(!exists(cdfname,where=where,inherits=FALSE)){
                      path <- .path.package("affy")
                      filename <- paste(cdfname,".rda",sep="")
                      load(file.path(path, "data", filename) ,
                           envir = where)
                    }
                    return(get(cdfname, envir=where))
                  }
                }
                
                if (what == "package") {
                  loc <- .find.package(cdfname, lib.loc=where, quiet=TRUE)
                  
                  if (identical(loc, character(0)))
                    next
                  ##stop(paste("AffyBatch: Looked for probes information in the package ", cdfname, "but could not find it.\n"))
                    
                  ##may be an option to try to autoload the package from
                  ##the bioconductor website woud be nice here
                  ## there is already a "how[[i]]$autoload" available
                  
                  if (length(loc) > 1)
                    warning(paste("several packages with a matching name. Using the one at", loc[1]))
                  
                  do.call("library", list(cdfname, lib.loc=dirname(loc[1])))
                  return(get(cdfname, envir=as.environment(paste("package:", cdfname, sep=""))))
                  ##object@cdfInfo <<- get(name, envir=as.environment(paste("package:",name, sep="")))
                }
              
                
                if (what == "file") {
                  ##now this is an actual Affymetrix filename
                  cdfname <- paste(object@cdfName,".CDF",sep="")
                  cdf <- read.cdffile(file.path(path.expand(where), cdfname))
                  ## ---> extra paranoia <---
                  if (cdf@cdfName != object@cdfName)
                  warning(paste("The CDF file identifies as", cdf@cdfName,
                                "while you probably want", object@cdfName))
                  ## ---> end <---
                  return(getLocations.Cdf(cdf))
                  ##object@cdfInfo <<- getLocations.Cdf(cdf)
                  rm(cdf)
                  gc() # since cdf can be rather large
                }
                
                if (what == "environment") {
                  if(exists(cleancdfname(object@cdfName),inherits=FALSE,where=where))
                  return(as.environment(get(cleancdfname(object@cdfName),
                                            where)))
                  ##object@cdfInfo <<- as.environment(get(name, where))
                }
              }
              stop(paste("AffyBatch: information about probe locations for ", object@cdfName, " could not be found"))
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

  
  ##show method
  if (debug.affy123) cat("--->show\n")
  setMethod("show", "AffyBatch",
            function(object) {
              
              cat("AffyBatch object\n")
              cat("size of arrays=", object@nrow, "x", object@ncol,
                  " features (", object.size(object), " Mb)\n", sep="")
              
              ## Location from cdf env
              try( cdf.env <- getCdfInfo(object) )
              if (! inherits(cdf.env, "try-error")) {
                num.ids <- length(ls(env=cdf.env))
              } else {
                warning("missing cdf environment !")
                num.ids <- "???"
              }
              
              cat("cdf=", object@cdfName,
                  " (", num.ids, " affyids)\n",
                  sep="")
              cat("number of experiments=",length(object),"\n",sep="")
              cat("number of genes=", length(geneNames(object)), "\n",sep="")
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
              
              i.probes <- match(which, c("pm", "mm"))
              ## i.probes will know if "[,1]" or "[,2]"
              
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
            function(object, genenames=NULL,xy=FALSE) 
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
            function(object, which=c("pm", "mm"),
                     genenames=NULL, LISTRUE=FALSE){
              
              which <- match.arg(which)
              
              index <- indexProbes(object, which, genenames)
              
              if(LISTRUE)
                ans <- lapply(index, function(i) object@exprs[i, ])
              else{
                index <- unlist(index)
                ans <- object@exprs[index, ]
                colnames(ans) <- sampleNames(object)
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
            function(object, genenames=NULL, LISTRUE=FALSE)
            probes(object, "pm", genenames, LISTRUE=LISTRUE),
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
  ##[[
  setMethod("[[", "AffyBatch",
            function(x, i, j, ...) { ##no need for j
              return(new("Cel",
                         intensity=matrix(intensity(x)[,i],x@ncol,x@nrow),
                         name=sampleNames(x)[i],
                         cdfName=x@cdfName))
              ##later we can get history from MIAME
            },where=where)
  
  ##[[ we need replacement that takes an entry by the Cel in value
    
  ##[ subseting. can only happen by sample. for now not by gene
  setMethod("[", "AffyBatch", function(x, i, j,..., drop=FALSE) {
    if( !missing(i) ) {
      phenoData(x) <- phenoData(x)[i, , ..., drop=FALSE]
      intensity(x) <- intensity(x)[ ,i, ..., drop=FALSE]
    }
    x
  },where=where)
  
  setReplaceMethod("[", "AffyBatch", function(x, i, j,..., value) {
    phenoData(x)[i,, ...] <- phenoData(value)[i, , ..., drop=FALSE]
    intensity(x)[,i]      <- intensity(value)[ ,i,... , drop=FALSE]
    x
  },where=where)

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
              object <- do.call(method, alist(object, ...))
              return(object)
            },
            where=where)

  
  ## expression value computation
  if (debug.affy123) cat("--->computeExprSet\n")
  if( !isGeneric("computeExprSet") )
    setGeneric("computeExprSet", function(x, bg.method, summary.method, ...)
               standardGeneric("computeExprSet"),
               where=where)
  
  setMethod("computeExprSet", signature(x="AffyBatch", bg.method="character", summary.method="character"),
            function(x, bg.method, summary.method, ids=NULL, verbose=TRUE, bg.param=list(), summary.param=list(), warnings=TRUE) {
              
              bg.method <- match.arg(bg.method, bg.correct.methods)
              summary.method <- match.arg(summary.method, express.summary.stat.methods)
              
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
              
              ##DEBUG: hackish (put global adjsutment names below
              if (bg.method %in% c("bg.correct.rma")) {
                if (verbose) cat("computing parameters for global background adjustement.....")
                all.l.pm.mat <- unlist(lapply(multiget(ids, env=getCdfInfo(x)),  function(x) if (ncol(x) == 2) x[,1]))
                all.l.mm.mat <- unlist(lapply(multiget(ids, env=getCdfInfo(x)),  function(x) if (ncol(x) == 2) x[,2]))
                all.param <- lapply(seq(1:n), function(i) {
                  notNA <- !(is.na(intensity(x)[, i][all.l.pm.mat]) | is.na(intensity(x)[, i][all.l.mm.mat]))
                  bg.parameters(intensity(x)[, i][notNA], intensity(x)[, i][notNA])
                })
                bg.param$all.param <- all.param
                if (verbose) cat(".....done.\n")
              }
              
             if (verbose) {
                cat(m, "ids to be processed\n")
                countprogress <- 0
              }
              
              ## loop over the ids
              mycall <- as.call(c(getMethod("express.summary.stat", signature=c("PPSet.container","character","character")),
                                  list(c.pps, method=summary.method, bg.correct=bg.method, param.bg.correct=bg.param, param.method=summary.param)))

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
                if (ncol(loc) == 2)
                  l.mm <- loc[ ,2]
                else
                  l.mm <- NA
                
                ## fill the PPSet.container
                ##c.pps@pmProbes <- matrix(NA, nrow=length(l.pm), ncol=n)
                ##c.pps@mmProbes <- matrix(NA, nrow=length(l.pm), ncol=n)

                np <- length(l.pm)
                
                ##names are skipped

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
                } else if (warnings) {
                  warning(paste("Error with affyid:", name.levels(cdf)[id]))
                }
                ## no need for an 'else' branching since exp.mat was initialized with NA
                
              }
              
              options(show.error.messages = TRUE)
              on.exit(NULL)
              
              if (verbose) cat("\n")

              ## instance exprSet
              ##if (verbose) cat("instancianting an exprSet.....")
              dimnames(exp.mat) <- list(ids, sampleNames(x))
              eset <- new("exprSet",
                          exprs=exp.mat,
                          se.exprs=matrix(), ##these needs to change
                          phenoData=phenoData(x),
                          description=description(x),
                          annotation=annotation(x),
                          notes=x@notes)
              ##if (verbose) cat(".....done.\n")
              
              return(eset)
            },
            where=where)
}

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
  







          

##this used to be in [[ and [[<-
##if we ever want to add other slots to this could
  ##be used in the "[[" method
  ##mysd <- sd(x)
  ##mymasks <- masks(x)
  ##myoutliers <- outliers(x)
  ##if statement make sure sd is actually somethings..
  ##same for other slots
  ##if(dim(mysd)[1]==ncols*nrows & dim(mysd)>=i)
  ##  mysd <- matrix(mysd[,i],ncols,nrows)
  ##if(dim(mymasks)[1]==ncols*nrows & dim(mymasks)>=i)
  ##  mymasks <- matrix(mymasks[,i],ncols,nrows)
  ##if(dim(myoutliers)[1]==ncols*nrows & dim(myoutliers)>=i)
  ##  myoutliers <- matrix(myoutliers[,i],ncols,nrows)
  
  ##cel <- new("Cel",
  ##           intensity=matrix(intensity(x)[i],ncols,nrows),
  ##           sd=mysd,
  ##           masks=mymaks,
  ##          outliers=myoutliers,
  ##           name=sampleNames(x)[i],
  ##           cdfName=x@cdfName,
  ##           history=history(x)[[i]])
  ##},
  ##where=where)
  
  ##this used to be part of "[[ 
  ##DEBUG: NA ?! watch this in next versions of R
              ## spotsd stuff to be really removed ?
              ##if (is.na.spotsd(x)) {
              #mysd <- matrix()
              ##} else {
              ##  mysd <- spotsd(x)[, , i]
              ##}
              ##oldim <- dim(intensity(x))
              ##dim(intensity(x)) <- c(x@nrow, x@ncol, x@nexp)
              ##new("Cel", intensity=intensity(x)[, , i], sd=mysd, name=sampleNames(x)[i], cdfName=x@cdfName, outliers=outliers(x)[[i]], masks=masks(x)[[i]], history=history(x)[[i]]) ## commented out becuase no 'outliers' or 'masks'.
              #cel <- new("Cel", intensity=intensity(x)[, , i], sd=mysd, name=sampleNames(x)[i], cdfName=x@cdfName, outliers=matrix(), masks=matrix(), history=history(x)[[i]])
              #dim(intensity(x)) <- oldim
              #return(cel)
            #},
            #where=where)

  ##if we ever want to add sd, masks, etc.. slots we can use this
  ## in the [[<- replacement method
  ##mysd <- sd(x)
  ##mymasks <- masks(x)
  ##myoutliers <- outliers(x)
  ##if statement make sure sd is actually somethings..
  ##if it is we put what should be put same for other slots
  ##if(dim(mysd)[1]==ncols*nrows & dim(mysd)>=i)
  ##  sd(x)[,i] <- as.vector(sd(value))
  ##if(dim(mymasks)[1]==ncols*nrows & dim(mymasks)>=i)
  ##  masks(x)[,i] <- as.vector(masks(value))
  ##if(dim(myoutliers)[1]==ncols*nrows & dim(myoutliers)>=i)
  ##  outliers(x)[,i] <- as.vector(outliers(value))
  ##},where=where)
  
  ##oldim <- dim(intensity(x))
  ##dim(intensity(x)) <- c(x@nrow, x@ncol, x@nexp)
  ##intensity(x)[,i] <- intensity(value)
#### spotsd ?
####if ((! is.na.spotsd(x)) & (spotsd(value) != c(NA)))
####  spotsd(x)[, , i] <- spotsd(value)
  ##x@name[i] <- sampleNames(value)
  ##
  ##if (x@cdfName != value@cdfName)
  ## warning("cdfName mismatch !")
  
  ##outliers(x)[[i]] <- outliers(value)
  ##masks(x)[[i]] <- masks(value)
  ##history(x)[[i]] <- history(value)
  ##dim(intensity(x)) <- oldim
  ##return(x)
  ##},
  ##where=where)
  
