.initAffyBatch <- function(where){
  
  if (debug.affy123) cat("-->initAffyBatch\n")

  ## Inherits from Affybatch
  ## The accessor 'intensity' gets what is in the slot 'exprs'
  setClass("AffyBatch",
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

#######################################################
### accessors
#######################################################

  if (debug.affy123) cat("--->accessors\n")
  
  ##intensity
  setMethod("intensity", signature(object="AffyBatch"),
            function(object) exprs(object),
            where=where)
  
  
  setReplaceMethod("intensity", signature(object="AffyBatch"),
                   function(object, value){
                     exprs(object) <- value
                     colnames(exprs(object)) <- sampleNames(object)
                     return(object)
                   }, where=where)

  ##for now, there is no accessor for se.exprs. we could use this to store
  ##sd, but if no one uses it... why do it
  
  setMethod("length",signature(x="AffyBatch"),
            function(x) ncol(exprs(x)), ##RI: assumes matrices
            where=where)
  if( !isGeneric("ncol") )
    setGeneric("ncol",where=where)
  setMethod("ncol",signature(x="AffyBatch"),
            function(x) x@ncol, ##RI: assumes matrices
            where=where)
  if( !isGeneric("nrow") )
    setGeneric("nrow",where=where)
  setMethod("nrow",signature(x="AffyBatch"),
            function(x) x@nrow, ##RI: assumes matrices
            where=where)

  
#######################################################
### methods
#######################################################
  
  
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
              second.try <- FALSE

              if (debug.affy123)
                cat("Trying to get cdfenv for", cdfname, "\n")
              
              
              i <- 0
              while(i < length(how)) {
                i <- i+1
                
                if (debug.affy123)
                  cat(i, ":")
                
                what <- how[[i]]$what
                where <- how[[i]]$where

                if (debug.affy123) {
                  cat("what=", what, "where=")
                  print(where)
                }
                
                if (what == "data") {
                  ##if we can get it from data dir. otherwise load package
                  
                  if(cdfname %in% do.call("data", list(package=where))$results[, 3]) {
                    ##RI: package="affy" doesnt work it has to be package=affy
                    ##    fix if you can
                    ##LG: weird stuff with data... but I had a workaround...
                    
                    where.env <- pos.to.env(match(paste("package:", where, sep = ""), search()))
                    
                    ## check if the cdfenv is already loaded. If not load it *in* the environment
                    ## of the package (where.env)
                    if( ! exists(cdfname, where = where.env, inherits = FALSE)) {
                      path <- .path.package(where)
                      filename <- paste(cdfname, ".rda", sep="")
                      load(file.path(path, "data", filename) ,
                           envir = where.env)
                    }
                    cdfenv <- get(cdfname, envir=where.env)
                    return(cdfenv)
                  }
                }
                
                if (what == "package") {
                  loc <- .find.package(cdfname, lib.loc=where, quiet=TRUE)
                  
                  if (!second.try && identical(loc, character(0))) {
                    ## before jumping to the next option, check the possibility to
                    ## download the missing cdfenv pack
                    
                    if (how[[i]]$autoload) {
                      cat(paste("Environment",cdfname,"is not available.\n"))
                      cat("This environment contains needed probe location information.\n\n")
                      
                      cat(paste("We will try to download and install the",
                                cdfname,"package.\n\n"))
                      if (! "package:reposTools" %in% search()) {
                        on.exit(detach("package:reposTools"))
                      }
                      
                      if (! require(reposTools))
                        stop("The package reposTools is required to download environments.\nPlease download and install it.\n")
                      
                      reposEntry <- getReposEntry(how[[i]]$repository)
                      
                      if (is.null(how[[i]]$installdir))
                        install.packages2(cdfname, reposEntry)
                      else
                        install.packages2(cdfname, reposEntry, how[[i]]$installdir)
                      
                      ## rewind the iterator i and try again
                      i <- i-1
                      second.try <- TRUE
                    }
                    ## jump to next way to get the cdfenv
                    next
                  }
                  if (length(loc) > 1)
                    warning(paste("several packages with a matching name. Using the one at", loc[1]))
                  
                  existsnow<- .find.package(cdfname, lib.loc=where, quiet=TRUE)

                  if(!identical(loc, character(0))){
                    do.call("library", list(cdfname, lib.loc=dirname(loc[1])))
                    
                    return(get(cdfname, envir=as.environment(paste("package:", cdfname, sep=""))))
                  }
                }
                
                
                if (what == "file") {
                  ##now this is an actual Affymetrix filename
                  cdfname <- paste(object@cdfName,".CDF",sep="")
                  cdf <- read.cdffile(file.path(path.expand(where), cdfname))
                  ## ---> extra paranoia <---
                  if (cdf@cdfName != object@cdfName)
                    warning(paste("The CDF file identifies as", cdf@cdfName,
                                  "while you probably want", object@cdfName))
                  ## ---> end of paranoia <---
                  return(getLocations.Cdf(cdf))
                }
                
                if (what == "environment") {
                  if(exists(object@cdfName,inherits=FALSE,where=where))
                    return(as.environment(get(object@cdfName,inherits=FALSE,envir=where)))
                }
              }
              stop(paste("\nWe could not find and/or install the necessary probe location information.\n",
                         "Here is a list of common problems and possible solutions:\n\n",
                         "Problem 1: You are not connected to the Internet.\n",
                         "Solution:  Try again once you connect to the Internet.\n\n",
                         "Problem 2: You do not have the necessary permissions to install packages.\n",
                         "Solution:  Ask your system administrator to install ",cleancdfname(object@cdfName), " package from:\n",
                         "           http://www.bioconductor.org/data/cdfenvs/cdfenvs.html\n\n",
                         "Problem 3: Necessary package not available from Bioconductor.\n",
                         "Solution:  Use makecdfenv package to create environment from CDF file.\n",
                         "           See the dealing_with_cdfenvs vignette for more details\n\n",
                         "NOTE: Once you install ",cleancdfname(object@cdfName)," you should not have this problem again.\n",
                         sep=""))
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

              ## Location from cdf env
              try( cdf.env <- getCdfInfo(object) )
              if (! inherits(cdf.env, "try-error")) {
                num.ids <- length(ls(env=cdf.env))
              } else {
                warning("missing cdf environment !")
                num.ids <- "???"
              }
              
              cat("AffyBatch object\n")
              cat("size of arrays=", nrow(object), "x", ncol(object),
                  " features (", object.size(object) %/% 1024, " kb)\n", sep="")
              cat("cdf=", object@cdfName,
                  " (", num.ids, " affyids)\n",
                  sep="")
              cat("number of samples=",length(object),"\n",sep="")
              cat("number of genes=", length(geneNames(object)), "\n",sep="")
              cat("annotation=",object@annotation,"\n",sep="")
              cat("notes=",object@notes,"\n",sep="")
            },
            where=where)

  
  if ( ! isGeneric("index2xy")) {
    setGeneric("indexProbes", function(object, which, ...)
               standardGeneric("indexProbes"), where=where)
  }

  
  ## indexProbes
  if( !isGeneric("indexProbes") )
    setGeneric("indexProbes", function(object, which, ...)
               standardGeneric("indexProbes"), where=where)

  setMethod("indexProbes", signature("AffyBatch", which="character"),
            function(object, which=c("pm", "mm","both"),
                     genenames=NULL, xy=FALSE) {
              
              which <- match.arg(which)
              
              i.probes <- match(which, c("pm", "mm", "both"))
              ## i.probes will know if "[,1]" or "[,2]"
              ## if both then [,c(1,2)]
              if(i.probes==3) i.probes=c(1,2)
              
              envir <- getCdfInfo(object)
              
              if(is.null(genenames)) 
                genenames <- ls(envir )
              
              ## shorter code, using the features of multiget
              ## (eventually more readable too)
              ## note: genenames could be confusing (the same gene can be
              ## found in several affyid (ex: the 3' and 5' controls)
              
              ans <-  multiget(genenames, pos, envir, iffail=NA)

              ## this kind of thing could be included in 'multiget' as
              ## an extra feature. A function could be specified to
              ## process what is 'multi'-get on the fly
              for (i in seq(along=ans)) {
                
                if ( is.na(ans[[i]][1]) )
                  next
                
                ##as.vector cause it might be a matrix if both
                tmp <- as.vector(ans[[i]][, i.probes])
                
                
                if (xy) {
                  warning("flag 'xy' is deprecated (because confusing)")
                  x <- tmp %% nrow(object)
                  x[x == 0] <- nrow(object)
                  y <- tmp %/% nrow(object) + 1
                  tmp <- cbind(x, y)
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
            function(object, genenames=NULL, xy=FALSE) 
            indexProbes(object, "pm", genenames=genenames, xy=xy),
            where=where
            )
  
  ##mmindex method
  if( !isGeneric("mmindex") )
    setGeneric("mmindex", function(object,...)
               standardGeneric("mmindex"), where=where)
  
  ##wrapper
  setMethod("mmindex", "AffyBatch",
            function(object,genenames=NULL, xy=FALSE) 
            indexProbes(object, "mm", genenames=genenames, xy=xy),
            where=where
            )                        

  
  ##probeNames method
  if( !isGeneric("probeNames") )
    setGeneric("probeNames", function(object, ...)
               standardGeneric("probeNames"), where=where)
  
  setMethod("probeNames","AffyBatch",
            function(object, genenames=NULL, mm=FALSE){
              if(mm) Index <- mmindex(object,genenames)
              else Index <- pmindex(object,genenames)
              reps <- unlist(lapply(Index,length))
              rep(names(Index),reps)
            },where=where)


  if( !isGeneric("probes") )
    setGeneric("probes", function(object, ...)
               standardGeneric("probes"), where=where)
  
  setMethod("probes", signature("AffyBatch"),
            function(object, which=c("pm", "mm"),
                     genenames=NULL, LISTRUE=FALSE, drop=FALSE){
              
              which <- match.arg(which)
              
              index <- indexProbes(object, which, genenames)
              
              if(LISTRUE)
                ans <- lapply(index, function(i) exprs(object)[i, ,drop=drop])
              else{
                index <- unlist(index)
                ans <- exprs(object)[index, ,drop=drop]
                colnames(ans) <- sampleNames(object)
                rownames(ans) <- names(index)
              }
              
              return(ans)
            },
            where=where)
  
  ##pm method
  if( !isGeneric("pm") )
    setGeneric("pm", function(object, ...)
               standardGeneric("pm"), where=where)
  
  setMethod("pm","AffyBatch",
            function(object, genenames=NULL, LISTRUE=FALSE)
            probes(object, "pm", genenames, LISTRUE=LISTRUE),
            where=where
            )
  
  if( !isGeneric("pm<-") )
    setGeneric("pm<-", function(object, value)
               standardGeneric("pm<-"), where=where)
  setReplaceMethod("pm", "AffyBatch", function(object, value){
    Dimnames <- dimnames(intensity(object))
    pmIndex <- unlist(pmindex(object))
    intensity(object)[pmIndex,] <- value
    dimnames(intensity(object)) <- Dimnames
    object
  }, where=where)
  


  ##mm method
  if( !isGeneric("mm") )
    setGeneric("mm", function(object, ...)
               standardGeneric("mm"), where=where)
  
  setMethod("mm",signature("AffyBatch"),
            function(object, genenames=NULL, LISTRUE=FALSE) probes(object, "mm", genenames, LISTRUE=LISTRUE),
            where=where
            )

  if( !isGeneric("mm<-") )
    setGeneric("mm<-", function(object, value)
               standardGeneric("mm<-"), where=where)
  setReplaceMethod("mm", "AffyBatch", function(object, value){
    Dimnames <- dimnames(intensity(object))
    mmIndex <- unlist(mmindex(object))
    intensity(object)[mmIndex,] <- value
    dimnames(intensity(object)) <- Dimnames
    object
  }, where=where)

###probeset
  if( !isGeneric("probeset") )
    setGeneric("probeset", function(object, ...)
               standardGeneric("probeset"), where=where)
  
  setMethod("probeset", "AffyBatch", function(object, genenames=NULL,
                                              locations=NULL){
    oldoptions <- getOption("BioC")
    
    if(is.null(locations)) ##use info in cdf
      envir <- getCdfInfo(object)
    else{
      ##if the user gives a list of locations let them use that as enviromnet
      envir <- new.env()
      multiassign(names(locations), locations, envir)
      object@cdfName <- "envir"
      newoptions <- oldoptions
      newoptions$affy$probesloc[[1]]$what <- "environment" 
      newoptions$affy$probesloc[[1]]$where <- parent.env(envir)
      options("BioC"=newoptions)
    }
    if(is.null(genenames))
      genenames <- ls(envir)
    
    p.pps <- vector("list", length(genenames))
    names(p.pps) <- genenames
    
    for (i in seq(along=genenames)) {
      
      i.pm <- indexProbes(object, "pm", genenames[i])[[1]]
      if (is.na(i.pm)[1])
        intensity.pm <- matrix()
      else
        intensity.pm <- intensity(object)[i.pm, , drop=FALSE]
      
      i.mm <- indexProbes(object, "mm", genenames[i])[[1]]
      if (is.na(i.mm)[1])
        intensity.mm <- matrix()
      else
        intensity.mm <- intensity(object)[i.mm, , drop=FALSE]
      
      p.pps[[i]] <- new("ProbeSet", id = genenames[i], pm = intensity.pm, mm = intensity.mm)
    }
    
    options("BioC"=oldoptions)
    return(p.pps)
  },where=where)
  
  if (debug.affy123) cat("--->[[\n")


  ##[[
  setMethod("[[", "AffyBatch",
            function(x, i, j, ...) { ##no need for j
              return(new("Cel",
                         intensity = matrix(intensity(x)[, i], ncol(x), nrow(x)),
                         name = sampleNames(x)[i],
                         cdfName = x@cdfName,
                         history = description(x)@preprocessing))
            },where=where)
  
  ##[[ we need replacement that takes an entry by the Cel in value
  
  ##[ subseting. can only happen by sample. for now not by gene
  setMethod("[", "AffyBatch", function(x, i, j,..., drop=FALSE) {
    if( !missing(i) ) {
      phenoData(x) <- phenoData(x)[i, , ..., drop=FALSE]
      intensity(x) <- intensity(x)[ ,i, ..., drop=FALSE]
    }
    return(x)
  },where=where)
  
  setReplaceMethod("[", "AffyBatch", function(x, i, j,..., value) {
    phenoData(x)[i,, ...] <- phenoData(value)[i, , ..., drop=FALSE]
    intensity(x)[,i]      <- intensity(value)[ ,i,... , drop=FALSE]
    return(x)
  },where=where)


  ## --- bg.correct

  if (debug.affy123) cat("--->bg.correct\n")

  if( !isGeneric("bg.correct") )
    setGeneric("bg.correct", function(object, method, ...)
               standardGeneric("bg.correct"), where=where)
  
  setMethod("bg.correct", signature(object="AffyBatch", method="character"),
            function(object, method=getOption("BioC")$affy$bgcorrect.method, ...) {

              ## simple for system to let one add background correction methods
              ## relies on naming convention
              
              method <- match.arg(method, bgcorrect.methods)
              
              methodname <- paste("bg.correct.", method, sep="")
              
              if (! exists(methodname))
                stop(paste("Unknown method (cannot find function", methodname, ")"))
              
              r <- do.call(methodname, alist(object, ...))
              
              return(r)
            }, where=where)


  ## --- normalize.methods
  if( !isGeneric("normalize.methods") )
    setGeneric("normalize.methods", function(object)
               standardGeneric("normalize.methods"),
               where=where)
  

  setMethod("normalize.methods", signature(object="AffyBatch"),
            function(object) {
              normalize.AffyBatch.methods
            },
            where=where)
  
  ## ---normalize  
  if (! isGeneric("normalize"))
    setGeneric("normalize", function(object, ...) standardGeneric("normalize"),
               where=where)
  
  setMethod("normalize", signature(object="AffyBatch"),
            function(object, method=getOption("BioC")$affy$normalize.method, ...) {
              method <- match.arg(method, normalize.AffyBatch.methods)
              if (is.na(method))
                stop("unknown method")
              method <- paste("normalize.AffyBatch", method, sep=".")
              object <- do.call(method, alist(object, ...))
              ## collect info in the attribute "normalization" 
              preproc <- c(description(object)@preprocessing,
                           list(normalization = attr(object, "normalization")))
              attr(object, "normalization") <- NULL
              ## and store it in MIAME
              MIAME <- description(object)
              MIAME@preprocessing <- preproc
              description(object) <- MIAME
              ##
              return(object)
            },
            where=where)

  
  ## --- expression value computation
  if (debug.affy123) cat("--->computeExprSet\n")
  if( !isGeneric("computeExprSet") )
    setGeneric("computeExprSet",
               function(x, pmcorrect.method, summary.method, ...)
               standardGeneric("computeExprSet"),
               where=where)
  
  setMethod("computeExprSet", signature(x="AffyBatch", pmcorrect.method="character", summary.method="character"),
            function(x, pmcorrect.method, summary.method, ids=NULL,
                     verbose=TRUE, summary.param=list(),
                     pmcorrect.param=list())
            {
              
              pmcorrect.method<- match.arg(pmcorrect.method, pmcorrect.methods)
              summary.method <- match.arg(summary.method, express.summary.stat.methods)
              
              n <- length(x)
              
              ## if NULL compute for all
              if (is.null(ids))
                ids <- geneNames(x)
              
              m <- length(ids)
              pps.warnings <- vector("list", length=m)
              
              ## cheap trick to (try to) save time
              c.pps <- new("ProbeSet",
                           pm=matrix(),
                           mm=matrix())
              
              ## matrix to hold expression values
              exp.mat <- matrix(NA, m, n)
              se.mat <- matrix(NA, m, n)
              
              if (verbose) {
                cat(m, "ids to be processed\n")
                countprogress <- 0
              }
              
              ## loop over the ids
              mycall <- as.call(c(getMethod("express.summary.stat",
                                            signature=c("ProbeSet","character", "character")),
                                  list(c.pps, pmcorrect=pmcorrect.method, summary=summary.method,
                                       summary.param=summary.param, pmcorrect.param=pmcorrect.param))
                                )
              ##only one character cause no more bg correct
              ##bg.correct=bg.method, param.bg.correct=bg.param,
              
              ##WHy not show error? took it out cause sometimes we
              ##get errors and couldnt see them.
              ##options(show.error.messages = FALSE)
              ##on.exit(options(show.error.messages = TRUE))
              
              CDFINFO <- getCdfInfo(x) ##do it once!
              
              for (i in seq(along=ids)) {
                
                id <- ids[i]
                
                if (verbose) {
                  if ( round(m/10) == countprogress) {
                    cat(".")
                    countprogress <- 0
                  }
                  else
                    countprogress <- countprogress + 1
                }
                ## locations for an id
                loc <- get(id, envir=CDFINFO)
                l.pm <- loc[, 1]
                if (ncol(loc) == 2)
                  l.mm <- loc[ ,2]
                else
                  l.mm <- integer()
                
                np <- length(l.pm)
                
                ##names are skipped

                c.pps@pm <- intensity(x)[l.pm, , drop=FALSE]
                c.pps@mm <- intensity(x)[l.mm, , drop=FALSE]
                
                ## generate expression values
                ## (wrapped in a sort of try/catch)
                mycall[[2]] <- c.pps
                ev <- try(eval(mycall))
                
                if (! inherits(ev, "try-error")) {
                  exp.mat[i, ] <- ev$exprs
                  se.mat[i,] <- ev$se.exprs
                  ## 
                } else {
                  pps.warnings[[i]] <- "Error"
                  ##warning(paste("Error with affyid:", id))
                }
                
              }
              
              #options(show.error.messages = TRUE)
              #on.exit(NULL)
              
              if (verbose) cat("\n")

              ## instance exprSet
              ##if (verbose) cat("instancianting an exprSet.....")
              dimnames(exp.mat) <- list(ids, sampleNames(x))
              dimnames(se.mat) <- list(ids, sampleNames(x))
              eset <- new("exprSet",
                          exprs=exp.mat,
                          se.exprs=se.mat,
                          phenoData=phenoData(x),
                          description=description(x),
                          annotation=annotation(x),
                          notes=c(notes(x)))
              ##if (verbose) cat(".....done.\n")

              attr(eset, "pps.warnings") <- pps.warnings
              return(eset)
              ##return(list(exprSet=eset, pps.warnings=pps.warnings))
            },
            where=where)


  ##some methods i was asked to add

  if( !isGeneric("image") )
    setGeneric("image",where=where)
  
  setMethod("image",signature(x="AffyBatch"),
            function(x, transfo=log, ...){
              scn <- prod(par("mfrow"))
              ask <- dev.interactive()
              which.plot <- 0
              for(i in 1:length(sampleNames(x))){
                which.plot <- which.plot+1;
                if(trunc((which.plot-1)/scn)==(which.plot-1)/scn && which.plot>1 && ask)  par(ask=TRUE)
                image(x[[i]], transfo=transfo, ...)
                par(ask=FALSE)}
            },where=where)
  

###boxplot
  if( !isGeneric("boxplot") )
    setGeneric("boxplot",where=where)
  setMethod("boxplot",signature(x="AffyBatch"),
            function(x,which="both",...){
              tmp <- description(x)
              if(class(tmp)=="MIAME") main <- tmp@title

              tmp <- unlist(indexProbes(x,which))
              tmp <- tmp[seq(1,length(tmp),len=5000)]

              boxplot(data.frame(log2(intensity(x)[tmp,])),main=main,range=0, ...)
            },where=where)

###hist
  if (debug.affy123) cat("--->hist\n")

  if( !isGeneric("hist") )
    setGeneric("hist",where=where)
  setMethod("hist",signature(x="AffyBatch"),function(x,...) plotDensity.AffyBatch(x,...),where=where)
  
}


##like for exprSet

"$.AffyBatch" <- function(affybatch, val)
    (pData(affybatch))[[as.character(val)]]

