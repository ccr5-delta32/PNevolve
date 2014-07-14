setwd("/run/user/1000/gvfs/smb-share:server=home.mpipz.mpg.de,share=homes,user=pieper/MPIPZ/cluster/PN_evolve")

stats <- read.table("PN-evo_MWU+KS_200_4.5.stats.tbl", header=F)

nstat <- subset(stats, stats[,2] == 0)
sstat <- subset(stats, stats[,2] == 1)

lnnMWU <- apply(nstat[,3:8], 2, function(x) log(x))
lnsMWU <- apply(sstat[,3:8], 2, function(x) log(x))

FishChi_n <- (apply(lnnMWU, 2, sum))
FishChi_n <- FishChi_n * -2
FishChi_np <- apply(FishChi_n, 1:2, function(x) pchisq(x, length(lnnMWU[,1])))
pchisq(FishChi_n[2], 2*length(lnnMWU[,1]),lower.tail=F)
test <- pchisq(136587.919, 18000,lower.tail=F)
length(lnnMWU[,1])

for (x in 3:8) {
	print(mean(nstat[,x]))
	print(mean(sstat[,x]))
}

Fisher.test <- function(p) {
  Xsq <- -2 * sum(log(p))
  p.val <- 1 - pchisq(Xsq, df = 2 * length(p))
  return(c(Xsq = Xsq, p.value = p.val))
} 

apply(sstat[,3:8], 2, Fisher.test)
apply(nstat[,3:8], 2, Fisher.test)


empidat <- read.csv("/run/user/1000/gvfs/smb-share:server=home.mpipz.mpg.de,share=homes,user=pieper/MPIPZ/x-perimentz/QTL_analysis/Meta/figures/acchisto/accpntn.csv", 
  header = T)
simdat  <- read.table("PN-evo_MWU+KS_200_4.5.resultsz.tbl", header=F)
nsimvec <- c(subset(simdat[,13:57], simdat[,2] == 0), recursive=T) 
ssimvec <- c(subset(simdat[,13:57], simdat[,2] == 1), recursive=T) 
empivec <- rep(empidat[,2], length(simdat[,1])/2)

ks.test(empidat[,2], nsimvec, alternative="two.sided")
ks.test(empidat[,2], nsimvec, alternative="greater")
ks.test(empidat[,2], nsimvec, alternative="less")
ks.test(empidat[,2], ssimvec, alternative="two.sided")
ks.test(empidat[,2], ssimvec, alternative="greater")
ks.test(empidat[,2], ssimvec, alternative="less")
ks.test(nsimvec, ssimvec, alternative="greater")
ks.test(nsimvec, ssimvec, alternative="less")

postscript(file = "ecdf_empi_sim.ps", width = 7.5, height = 7.5, paper = "special", 
  horizontal = F)
plot.ecdf(empidat[, 2], cex = 0.4, xlim = c(0, 4.05), main = "Cummulative distribution functions", 
  xlab = "mean petal number", ylab = "fraction of observations", las = 1)
par(new = T)
plot.ecdf(ssimvec, col = 4, cex = 0.4, xlim = c(0, 4.05), axes=F, main=NULL, xlab="", ylab="")
par(new = T)
plot.ecdf(nsimvec, col = 2, cex = 0.4, xlim = c(0, 4.05), axes=F, main=NULL, xlab="", ylab="")
legend(0,0.98, legend=c("empirical data","simulated data under neutrality","simulated data under selection"), col=c(1,2,4), fill=c(1,2,4))
dev.off() 

?ks.test

length(nsimvec)
length(ssimvec)

postscript(file="simhis_neutral.ps", width=6, height=6, paper="special", horizontal=F)
par(mar=c(5.1,5.1,4.1,1.1), mgp=c(4,1,0))
options(scipen=5)
hist(nsimvec, breaks=seq(0,4,by=0.4), col="grey52", main="", xlab="mean petal number", ylab="Frequency (number of simulated genotypes)", las=1, cex.lab=1.2)
box()
dev.off()
postscript(file="simhis_selection.ps", width=6, height=6, paper="special", horizontal=F)
par(mar=c(5.1,5.1,4.1,1.1), mgp=c(4,1,0))
hist(ssimvec, breaks=seq(0,4,by=0.4), col="grey52", main="", xlab="mean petal number", ylab="Frequency (number of simulated genotypes)", las=1, cex.lab=1.2) 
box()
dev.off()
options()$scipen
empidat[,2]
head(simdat)

# Asis' mutationevent function 
mutationEvent <- function(n.qtls, base.prob=7.5e-5 ) { p.mut <- runif( 1, 0, 1 ); p.mut <= n.qtls * base.prob }
