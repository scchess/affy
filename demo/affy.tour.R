## A quick demo to overview the package
##         -- Laurent

if(dev.cur() <= 1) get(getOption("device"))()

opar <- par(ask= (interactive() &&
                  (.Device %in% c("X11", "GTK", "windows", "Macintosh"))),
            bg="cornsilk", mfrow=c(1,1))



## load the data

data(listcel)
data(CDF.example)


## display the image of the data in the CEL file

image(listcel[[1]],transfo=log)


image(listcel[[1]],transfo=log)

## find the locations for probes corresponding to a given ID

l.pm <- locate.name("AFFX-BioC-5_at", CDF.example, type="pm")
plot.location(l.pm, CDF.example, col="red", pch=16)
l.mm <- locate.name("AFFX-BioC-5_at", CDF.example, type="mm")
plot.location(l.mm, CDF.example, col="blue", pch=16)

#legend(0.4,0,c("perfect match","mismatch"),c("red","blue"),bg="white")

rm(l.pm, l.mm)



## The '5', 'M' or '3' in the names means the set of probes relates to
## the 5-prime, middle or 3-prime sectors  respectively.
## The at or st ending means the sequence relates to the complementary
## sequence of the gene or not respectively.

namesspot <- c("AFFX-BioC-5_at","AFFX-BioC-3_at")

p <- lapply(namesspot, get.PPSet, CDF.example, listcel[[1]])

ylimi <- range(unlist(lapply(p, function(x) x@probes)), 0)

plot.new()
par(mfrow=c(2,2))

for(i in 1:length(namesspot)) barplot(p[[i]],ylim=ylimi)

rm(ylimi,p)

plot(1:10,1:10,type="n", yaxt="n", xlab="", xaxt="n",ylab="")
legend(3,8,c("PM","MM"),fill=c("red","blue"))
text(3,c(4,3,2),
     c("Probe pair sets for","the 5' and 3' ends","of a control gene."))


## normalize
plot.new()
par(mfrow=c(2,2))

nat <- pmormm(CDF.example)

plot(listcel[[1]]@intensity, listcel[[2]]@intensity, xlab="CEL file 1", ylab="CEL file 2",main="raw values",sub="all probes plotted",type="n")
points(listcel[[1]]@intensity[nat], listcel[[2]]@intensity[nat], col="red")
points(listcel[[1]]@intensity[!nat], listcel[[2]]@intensity[!nat], col="blue")
points(listcel[[1]]@intensity[is.na(nat)], listcel[[2]]@intensity[is.na(nat)], pch="+")
legend(25000, 15000, c("PM","MM","Unknown","identity line"), c("red","blue","black","grey"), bg="white")
lim <- range(par()$usr)
points(lim,lim,type="l",col="gray")
rm(nat)

listcel.n <- normalize(listcel, f.cdf=CDF.example, method="constant", refindex=2)
plot(listcel.n[[1]]@intensity, listcel.n[[2]]@intensity, xlab="CEL file 1", ylab="CEL file 2",main="normalized by constant",sub="all probes plotted")
lim <- range(par()$usr)
points(lim,lim,type="l",col="gray")


listcel.n <- normalize(listcel, f.cdf=CDF.example, method="invariantset")
i.set <- listcel.n[[2]]@history$invariantset

plot(listcel[[1]]@intensity, listcel[[2]]@intensity, xlab="CEL file 1", ylab="CEL file 2",main="raw values",sub="all probes plotted")
points(listcel[[1]]@intensity[i.set], listcel[[2]]@intensity[i.set], col="orange",pch=16)
lim <- range(par()$usr)
points(lim,lim,type="l",col="gray")
points(smooth.spline(listcel[[1]]@intensity[i.set], listcel[[2]]@intensity[i.set]),type="l",col="red")
lim <- range(par()$usr)
points(lim,lim,type="l",col="gray")
legend(25000,15000,c("invariant set","identity line","spline through the invariant set"),c("orange","grey","red"),bg="white")

plot(listcel.n[[1]]@intensity, listcel.n[[2]]@intensity, xlab="CEL file 1", ylab="CEL file 2",main="normalized by invariant set",sub="all probes plotted")
points(listcel.n[[1]]@intensity[i.set], listcel.n[[2]]@intensity[i.set], col="orange",pch=16)
lim <- range(par()$usr)
points(lim,lim,type="l",col="gray")
legend(20000,10000,c("invariant set","identity line"),c("orange","grey"),bg="white")

##
## Normalization and its effects can be observed in a couple of commands
##

par(mfrow=c(2,2))
data(Dilution)
boxplot(Dilution)
boxplot(normalize(Dilution))


##
# generates pseudo repeated measurments
# by adding noise to two differents measurements
##
p <- new("PPSet.container", x=vector("list", length=4), content="PPSet")
p[[1]] <- get.PPSet("A28102_at",CDF.example,listcel[[1]])
p[[2]] <- get.PPSet("A28102_at",CDF.example,listcel[[2]])
p[[3]] <- p[[1]]
p[[4]] <- p[[2]]

p[[3]]@probes$pm[c(13,7)] <- p[[3]]@probes$pm[c(13,7)] + c(500,-400)
p[[4]]@probes$pm[c(3,7)] <- p[[4]]@probes$pm[c(3,7)] + runif(2,10,2000)

ylimi <- range(unlist(lapply(p@x, function(y) y@probes)), 0)

#par(mfcol=c(4,2))
nmet <- 2

layout(matrix(c(1:4,rep(5:(5+nmet-1), times=rep(4,nmet))),4,1+nmet), width=c(4,rep(1,nmet)))

for (i in 1:4) barplot(p[[i]],ylim=ylimi, main=paste(p[[i]]@name, " - hybridization ", i))

ev <- generateExprVal.PPSet.container(p,method="playerout")
barplot(rev(c(ev)),main="expression values using\nplayerout",
        names.arg=rev(paste("hybrid. ",1:4)), horiz=TRUE)
ev <- generateExprVal.PPSet.container(p,method="avgdiff")
barplot(rev(c(ev)),main="expression values using\navgdiff",
        names.arg=rev(paste("hybrid. ",1:4)), horiz=TRUE)
#ev <- generateExprVal(p,method="affy")
#barplot(c(ev),main="expression values using\nAffymetrix truncated average difference")
#ev <- generateExprVal(p,method="liwong.reduced")
#barplot(c(ev),main="expression values using\nLi and Wong reduced model")

par(opar)
rm(opar)



