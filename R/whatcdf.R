whatcdf <- function(filename, compress=getOption("BioC")$affy$compress.cel){
  
  ##finds what cdf environment to use with cdf file
  tmp <- getInfoInAffyFile(filename,"CEL","HEADER","DatHeader",compress=compress) ##find appropriate line
  tmp <- strsplit(tmp," ")[[1]] #split by space
  tmp <- tmp[grep(".1sq",tmp)] #pick the one with 1sq (from experience)
  tmp <- gsub("\.1sq","",tmp) #take out .1sq
  tmp <- gsub("_","",tmp) #take out underscore
  tmp <- tolower(tmp) #make lower case
  return(tmp)
}
