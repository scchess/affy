getCdfInfo <- function(object,  how=getOption("BioC")$affy$probesloc, verbose=TRUE) {
    ## cdfname is the cdf environment
    ## methods is a vector detailing how to get the file - one of
    ## 'environment', 'data', 'library', 'bioC'
    ## where is used for specific information corresponding to method

    if (length(how) == 0)
        stop("No available method to obtain CDF file")

    cdfname <- cleancdfname(cdfName(object))

    badOut <- list()
    for (i in 1:length(how)) {
        cur <- how[[i]]
        out <- switch(cur$what,
                      "environment" = cdfFromEnvironment(cdfname,
                      cur$where, verbose),
                      "data" = cdfFromData(cdfname, cur$where, verbose),
                      "libPath" = cdfFromLibPath(cdfname, cur$where,
                      verbose=TRUE),
                      "bioC" = cdfFromBioC(cdfname, cur$where,
                      verbose)
                      )

        if (is.environment(out))
            return(out)
        else
            badOut <- c(badOut, out)
    }

    stop(paste("Could not obtain CDF environment, problems encountered:",
               paste(unlist(badOut),collapse="\n"),sep="\n"))
}

cdfFromData <- function(cdfname, pkg, verbose=TRUE) {
    if (verbose)
        print(paste("Attempting to locate",cdfname,
                    "in the data directory of package", pkg))

    if(cdfname %in% do.call("data", list(package=pkg))$results[, 3]) {
        where.env <- pos.to.env(match(paste("package:", pkg, sep = ""), search()))

        ## check if the cdfenv is already loaded. If not load it *in* the environment
        ## of the package (where.env)
        if(!exists(cdfname, where = where.env, inherits = FALSE)) {
            path <- .path.package(pkg)
            filename <- paste(cdfname, ".rda", sep="")
            load(file.path(path, "data", filename) ,
                 envir = where.env)
        }
        cdfenv <- get(cdfname, envir=where.env)
        return(cdfenv)
    }
    else
        return(list(paste("Data for package",pkg,"did not contain",cdfname)))
}

cdfFromEnvironment <- function(cdfname, where, verbose=TRUE) {
    if (verbose)
        print(paste("Attempting to locate",cdfname,"in specified environment"))

    if (exists(cdfname, inherits=FALSE, where=where))
        return(as.environment(get(cdfname,inherits=FALSE,envir=where)))
    else {
        if (verbose)
            print(paste("Specified environment does not contain",cdfname))
        return(list(paste("Specified environment specified did not contain",cdfname)))
    }
}

cdfFromBioC <- function(cdfname, lib=.libPaths()[1], verbose=TRUE) {
    require(reposTools) || stop("Package 'reposTools' is required",
                                " for this operation.")

    if (verbose)
        print(paste("Attempting to obtain",cdfname,"from Bioconductor website"))

    ## First search the user's libPaths to see if it is installed
    if (verbose)
        print(paste("Checking to see if the environment",
                    cdfname,"is already installed ..."))
    if (is.installed(cdfname)) {
        if (verbose)
            print(paste("The environment",cdfname,"is already installed."))
    }
    else {
        if (verbose)
            print(paste("The environment ",cdfname," was not found in",
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
                    return(list("Bioconductor - lib does not exist"))
                }
            }

            if (file.access(lib,mode=2) < 0) {
                if (verbose) {
                    print(paste("You do not have write access to",lib,
                               "\nPlease check your permissions or provide",
                               "a different 'lib' parameter"))
                    return(list("Bioconductor - lib is not writeable"))
                }
            }

            z <- install.packages2(cdfname, lib=lib)
            if(! cdfname %in% updatedPkgs(z)) {
                if (verbose)
                    print(paste("Environment",cdfname,
                                "was not found in the Bioconductor",
                                "repository."))
                return(list(paste("Bioconductor -",cdfname,"not available")))
            }
            else
                if (verbose)
                    print(paste("Installation of environment",
                                cdfname, "was succesful."))
        }
        else {
            if (verbose)
                print(paste("The current operation could not access",
                            "the Bioconductor repository.  Please",
                            "check your internet connection, and",
                            "report further problems to",
                            "bioconductor@stat.math.ethz.ch"))
            return(list("Bioconductor - could not connect"))
        }
    }

    ## Now load the library and return the environment
    do.call("library", list(cdfname, lib.loc=lib))
    ## !!! Double check that library is actually loaded
    if (! cdfname %in% .packages()) {
        ## package was not properly loaded
        if (verbose)
            print(paste("The package", cdfname,
                        "could not be loaded"))
        return(list("Bioconductor - package downloaded but not loadable"))
    }
    else
        return(get(cdfname,
                   envir=as.environment(paste("package:", cdfname, sep=""))))
}

cdfFromLibPath <- function(cdfname, lib = NULL, verbose=TRUE) {
    ## First check to see if package is installed
    if (verbose)
        print(paste("Checking to see if package",cdfname,
                    "is already installed"))

    if (length(.find.package(cdfname, lib.loc=lib, quiet=TRUE)) == 0)
        return(list(paste("Library - package",cdfname,"not installed")))

    ## See if package is already loaded
    if (cdfname %in% .packages()) {
        if (verbose)
            print(paste("The package", cdfname, "is already loaded"))
    }
    else {
        if (verbose)
            print(paste("Attempting to load package", cdfname))
        ## Attempt to load the library requested
        do.call("library", list(cdfname, lib.loc=lib))

        ## Check to see if it got loaded
        if (! cdfname %in% .packages()) {
            ## package didn't get loaded
            if (verbose)
                print(paste("The package", cdfname, "could not be loaded"))
            return(list(paste("Library - package",cdfname,"is not loadable")))
        }
    }

    return(get(cdfname, envir=as.environment(paste("package:", cdfname, sep=""))))
}
