getCdfInfo <- function(object, how=getOption("BioC")$affy$probesloc) {
    ## "how" is a list. each element of the list must have an element
    ## tagged "what" and an element tagged "where"
    ## "what" can be "package" or "environment"
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
            ##need this until cdfenv.example is better handled
            ##if we can get it from data dir. otherwise load package

            if(cdfname %in% do.call("data", list(package=where))$results[, 3]) {
                where.env <- pos.to.env(match(paste("package:", where, sep = ""), search()))

                ## check if the cdfenv is already loaded. If not load it *in* the environment
                ## of the package (where.env)
                if(!exists(cdfname, where = where.env, inherits = FALSE)) {
                    path <- .path.package(where)
                    filename <- paste(cdfname, ".rda", sep="")
                    load(file.path(path, "data", filename) ,
                         envir = where.env)
                }
                cdfenv <- get(cdfname, envir=where.env)
                return(cdfenv)
            }
            next
        }

        if (what == "package") {
            loc <- .find.package(cdfname, lib.loc=where, quiet=TRUE)

            if (!second.try && identical(loc, character(0))) {
                ## before jumping to the next option, check the possibility to
                ## download the missing cdfenv pack

                if (how[[i]]$autoload) {
                    cat(paste("Environment",cdfname,"is not available.\n"))
                    cat("This environment contains needed probe location information.\n")

                    cat(paste("We will try to download and install the",
                              cdfname,"package.\n\n"))

                    if (is.null(how[[i]]$installdir))
                        status <- getCDFenv(cdfname,verbose=TRUE)
                    else
                        status <- getCDFenv(cdfname,
                                            how[[i]]$installdir,
                                            verbose=TRUE)

                    if (status) {
                        ## rewind the iterator i and try again
                        i <- i-1
                        second.try <- TRUE
                    }
                }
                ## jump to next way to get the cdfenv
                next
            }

            if (length(loc) > 1)
                warning(paste("several packages with a matching name. Using the one at", loc[1]))

            if(!identical(loc, character(0))){
                do.call("library", list(cdfname, lib.loc=dirname(loc[1])))
                return(get(cdfname, envir=as.environment(paste("package:", cdfname, sep=""))))
            }
            next
        }

        if (what == "environment") {
            if(exists(object@cdfName,inherits=FALSE,where=where))
                return(as.environment(get(object@cdfName,inherits=FALSE,envir=where)))
                next
        }
    }
}


getCDFenv <- function(env, lib=.libPaths()[1], verbose=TRUE) {
    require(reposTools) || stop("Package 'reposTools' is required",
                                " for this operation.")

    ## !! First search the user's libPaths to see if it is installed
    if (verbose)
        print(paste("Checking to see if the environment",
                    env,"is already installed ..."))
    if (is.installed(env)) {
        if (verbose)
            print(paste("The environment",env,"is already installed."))
    }
    else {
        if (verbose)
            print(paste("The environment ",env," was not found in",
                        " these directories: ",
                        paste(.libPaths(), collapse=", "),
                        ".  Now searching the internet repository.",
                        sep=""))
        if (verbose)
            print(paste("Checking to see if your internet",
                        "connection works ..."))
        if (testBioCConnection()) {

            ## Check for file permissions
            if (file.access(lib, mode=0) < 0) {
                if (verbose) {
                    print(paste("Directory",lib,"does not seem to exist.\n",
                                "Please check your 'lib' parameter and try again"))
                    return(invisible(FALSE))
                }
            }

            if (file.access(lib,mode=2) < 0) {
                if (verbose) {
                    print(paste("You do not have write access to",lib,
                               "\nPlease check your permissions or provide",
                               "a different 'lib' parameter"))
                    return(invisible(FALSE))
                }
            }

            z <- install.packages2(env, lib=lib)
            if(! env %in% updatedPkgs(z)) {
                if (verbose)
                    print(paste("Environment",env,
                                "was not found in the Bioconductor",
                                "repository."))
                return(invisible(FALSE))
            }
            else
                if (verbose)
                    print(paste("Installation of environment",
                                env, "was succesful."))
        }
        else {
            if (verbose)
                print(paste("The current operation could not access",
                            "the Bioconductor repository.  Please",
                            "check your internet connection, and",
                            "report further problems to",
                            "bioconductor@stat.math.ethz.ch"))
            return(invisible(FALSE))
        }
    }
    return(invisible(TRUE))
}

testBioCConnection <- function() {
    ## Stifle the "connected to www.... garbage output
    curNetOpt <- getOption("internet.info")
    on.exit(options(internet.info=curNetOpt), add=TRUE)
    options(internet.info=3)

    ## First check to make sure they have HTTP capability.  If they do
    ## not, there is no point to this exercise.
    http <- as.logical(capabilities(what="http/ftp"))
    if (http == FALSE)
        return(FALSE)

    ## find out where we think that bioC is
    bioCoption <- getOption("BIOC")
    if (is.null(bioCoption))
        bioCoption <- "http://www.bioconductor.org"

    ## Now check to see if we can connect to the BioC website
    biocURL <- url(paste(bioCoption,"/main.html",sep=""))
    options(show.error.messages=FALSE)
    test <- try(readLines(biocURL)[1])
    options(show.error.messages=TRUE)
    if (inherits(test,"try-error"))
        return(FALSE)
    else
        close(biocURL)

    return(TRUE)
}
