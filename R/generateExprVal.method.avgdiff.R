## Currently, the input is a list of data.frames.
## All the elements in the list correspond to the same gene
## Each data.frame represent a PPSet from an array, and has two
## columns one is called 'pm' and the other 'mm'

generateExprVal.method.avgdiff <- function(probes) {
  unlist(lapply(probes,function(x){mean(x$pm - x$mm)}))
}
