## A quick demo to overview the package
##         -- Laurent

if(dev.cur() <= 1) get(getOption("device"))()

opar <- par(ask= (interactive() &&
                  (.Device %in% c("X11", "GTK", "windows", "Macintosh"))),
            bg="cornsilk", mfrow=c(1,1))



## load the data

data(CDF.example)
data(affybatch.example)

## display the image of the data in the CEL file

cel <- affybatch.example[[1]]

image(cel,transfo=log)

image(cel,transfo=log)

## find the locations for probes corresponding to a given ID

l.pm <- locate.name("AFFX-BioC-5_at", CDF.example, type="pm")
plotLocation(l.pm, CDF.example, col="red", pch=16)
l.mm <- locate.name("AFFX-BioC-5_at", CDF.example, type="mm")
plotLocation(l.mm, CDF.example, col="blue", pch=16)

#legend(0.4,0,c("perfect match","mismatch"),c("red","blue"),bg="white")

rm(l.pm, l.mm)



## The '5', 'M' or '3' in the names means the set of probes relates to
## the 5-prime, middle or 3-prime sectors  respectively.
## The at or st ending means the sequence relates to the complementary
## sequence of the gene or not respectively.

namesspot <- c("AFFX-BioB-5_at","AFFX-BioB-M_at", "AFFX-BioB-3_at")

p <- probeset(affybatch.example, genenames=namesspot)

par(mfrow=c(3,3))

for (pps in p)
  barplot(pps)


## normalize
plot.new()
par(mfrow=c(2,2))

nat <- pmormm(CDF.example)

cel2 <- affybatch.example[[2]]
plot(intensity(cel), intensity(cel2), xlab="CEL file 1", ylab="CEL file 2",main="raw values",sub="all probes plotted",type="n")
points(intensity(cel)[nat], intensity(cel2)[nat], col="red")
points(intensity(cel)[!nat], intensity(cel2)[!nat], col="blue")
points(intensity(cel)[is.na(nat)], intensity(cel2)[is.na(nat)], pch="+")
legend(25000, 15000, c("PM","MM","Unknown","identity line"), c("red","blue","black","grey"), bg="white")
abline(0, 1, type="l", col="gray")
rm(nat)

abatch.n <- normalize(affybatch.example, method="constant", refindex=2)
plot(intensity(abatch.n[[1]]), intensity(abatch.n[[2]]), xlab="CEL file 1", ylab="CEL file 2",main="normalized by constant",sub="all probes plotted")
abline(0, 1, type="l", col="gray")


abatch.n <- normalize(affybatch.example, method="invariantset")
i.set <- history(abatch.n[[1]])$invariantset

plot(intensity(cel), intensity(cel2), xlab="CEL file 1", ylab="CEL file 2",main="raw values",sub="all probes plotted")

abline(0, 1, type="l", col="gray")
legend(25000,15000,c("invariant set","identity line","spline through the invariant set"),c("orange","grey","red"),bg="white")

plot(intensity(abatch.n[[1]]), intensity(abatch.n[[2]]), xlab="CEL file 1", ylab="CEL file 2",main="normalized by invariant set",sub="all probes plotted")
points(intensity(abatch.n[[1]])[i.set], intensity(abatch.n[[2]])[i.set], col="orange",pch=16)
abline(0, 1, type="l", col="gray")
legend(20000,10000,c("invariant set","identity line"),c("orange","grey"),bg="white")

##
## Normalization and its effects can be observed in a couple of commands
##

#par(mfrow=c(2,2))
#data(Dilution)
#boxplot(Dilution)
#boxplot(normalize(Dilution))


p <- probeset(affybatch.example, "A28102_at")

#par(mfcol=c(4,2))
#mymethods <- c("avgdiff", "playerout", "liwong")
#nmet <- length(mymethods)

#layout(matrix(c(1:4, rep(5:(5+nmet-1), times=rep(4, nmet))), 4, 1+nmet),
#       width=c(4, rep(1, nmet)))

#for (i in 1:4) barplot(p[[i]], ylim=ylimi, main=paste(p[[i]]@name, " - hybridization ", i))

#for (i in 1:nmet) {
#  ev <- express.summary.stat(p, method=mymethods[i], bg.correct="bg.correct.pmonly")
#  barplot(rev(c(ev)),main=paste("expression values using\n",mymethods[i], sep=""),
#          names.arg=rev(paste("hybrid. ",1:4)), horiz=TRUE)
#}
#par(opar)
#rm(opar)
  



