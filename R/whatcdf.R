##this function changes the affymetrix cdf file name to the Bioconductor
##annotation name for that cdf file
## note: we had a hard time finding exact rules to match what is in the
## CEL file with what is in the CDF file
## ex: CEL says 'ecoli' while CDF says 'ecoligenome'
## or: CEL says '' while CDF says hu6800.1sq
cleancdfname <- function(cdfname, addcdf=TRUE) {
  i <- match(cdfname, mapCdfName$inCDF)
  if (is.na(i)) {
    tmp <- tolower(cdfname) #make lower case
    tmp <- gsub("_", "", tmp) #take out underscore
    tmp <- gsub("-", "", tmp) #take out underscore
    tmp <- gsub("\ ", "", tmp) ##take out spaces
    if(addcdf) tmp <- paste(tmp, "cdf", sep="")
  } else {
    tmp <- mapCdfName$inBioC[1]
  }
  return(tmp)
}
##this funnction gets the cdf from a celfile

whatcdf <- function(filename, compress=getOption("BioC")$affy$compress.cel){
  
  ##finds what cdf environment to use with cdf file
  tmp <- getInfoInAffyFile(filename,"CEL","HEADER","DatHeader",compress=compress) ##find appropriate line
  tmp <- strsplit(tmp," ")[[1]] #split by space
  tmp <- tmp[grep(".1sq",tmp)] #pick the one with 1sq (from experience)
  if (identical(tmp, character(0))) {
    warning("could not find CDF name, setting it to 'unknown'")
    tmp <- "unknown"
  }
  else {
    tmp <- gsub("\.1sq","",tmp) #take out .1sq
  }
  return(tmp)
}

  
