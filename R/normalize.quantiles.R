.First.lib <- function(lib, pkg) library.dynam("affy", pkg, lib) 

normalize.quantiles <- function(x){

  rows <- dim(x)[1]
  cols <- dim(x)[2]
  
  matrix(.C("qnorm_c",as.double(as.vector(x)),as.integer(rows),as.integer(cols))[[1]],rows,cols)
}


