##this function changes the affymetrix cdf file name to the Bioconductor
##annotation name for that cdf file
cleancdfname <- function(cdfname,addcdf=TRUE){
  tmp <- gsub("_","",cdfname) #take out underscore
  tmp <- gsub("-","",tmp) #take out underscore
  tmp <- gsub("\ ","",tmp) ##take out spaces
  tmp <- tolower(tmp) #make lower case
  if(addcdf) tmp <- paste(tmp,"cdf",sep="")
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
