library(affy)

## CEL file

n <- 50

## --- ripped from the example in 'persp()'
x <- seq(-10, 10, length=n)
y <- x
f <- function(x,y)
  {
    r <- sqrt(x^2+y^2)
    10 * sin(r)/r
  }
z <- outer(x, y, f)
z[is.na(z)] <- 1
## ------------------
z <- z - min(z)
 
#DEBUG: cleaner way to get a temp file ?
tmpfile <- .getTmpFileName()
cel <- new("Cel", intensity=z, sd=z/10, name="", cdfName="dummy.1sq",
            masks=matrix(c(1,1), nrow=1), outliers=matrix(c(1,1), nrow=1))

write.celfile(cel, tmpfile)

cat("---> read.affybatch...\n")
afbatch <- read.affybatch(filenames=c(tmpfile, tmpfile))
cat("done.\n")

unlink(tmpfile)


## fake environment

dummy <- new.env(hash=T)
index <- cbind(runif(10, 1, n), runif(10, 1, n))
assign("gene.a", index, envir=dummy)
index <- cbind(runif(10, 1, n), runif(10, 1, n))
assign("gene.b", index, envir=dummy)

## set options to look for the locations of the probes in the
## current environment
opt <- getOption("BioC")
old.opt <- opt
opt$affy$probesloc.what="environment"
opt$affy$probesloc.where=environment()
options(BioC=opt)

## get a Cel in the AffyBatch
cat("---> getting a Cel from an AffyBatch...\n")
cel <- afbatch[[1]]
cat("done.\n")

## normalize the AffyBatch
cat("---> normalizing an AffyBatch...\n")
n.afbatch <- normalize(afbatch, method="constant")
cat("done.\n")

## compute expression values
##cat("---> normalizing an AffyBatch...\n")
#eset <- generateExprSet(n.afbatch, method="liwong")
##cat("done.\n")


## restore old options
options(BioC=old.opt)
