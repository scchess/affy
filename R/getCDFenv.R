getCDFenv <- function(env, verbose=TRUE) {
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
            z <- install.packages2(env)
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
