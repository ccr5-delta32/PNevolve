# Bjorn Pieper. MPIPZ, Cologne, Germany. Make pretty graphs of s1VSs2 scan results

setwd("/home/bpuser/server/home/cluster/PN_evolve/run2_")

scanmu <- read.table("run2_v1_s2VSs1.mu.tbl", nrows=256, header=F)

xs1 <- seq(min(scanmu[,2]), max(scanmu[,2]), by=1)
xs2 <- seq(min(scanmu[,3]), max(scanmu[,3]), by=0.01)

indices <- cbind(seq(max(scanmu[,3]), min(scanmu[,3]), by=-0.1), seq(max(scanmu[,2]), min(scanmu[,2]), by=-10))

magicnumber <- 0.85
png(file="scan_s1VSs2_results.png", width=3500, height=2000, unit="px", res=300)

for (x in 1:length(indices[,1])) {
  if (x == 1) {
    ylim <- c(0,(1.01*max(scanmu[,4]+max(scanmu[,5]))))
    red  <- 255; green <- 0; blue <- 0
    colo <- rgb(red, green,blue,255,maxColorValue=255)
  } else {
    red = red - (255/length(seq(max(scanmu[,3]), min(scanmu[,3]), by=-0.1)))
    green = green + (255/length(seq(max(scanmu[,3]), min(scanmu[,3]), by=-0.1)))
    blue = blue + 0.7* (255/length(seq(max(scanmu[,3]), min(scanmu[,3]), by=-0.1)))
    colo <- rgb(red, green,blue,255,maxColorValue=255)
  }
  s1temp <- subset(scanmu, scanmu[,3] == as.character(indices[x,1]) )
  s1res  <- s1temp[,4]
  S1   <- s1temp[,2]
  s1fit<- nls(s1res ~ (A / S1) + (B + log(S1, C)), start=list(A=-5, B=1.3, C=100))
  s2temp <- subset(scanmu, scanmu[,2] == as.character(indices[x,2]))
  s2res  <- s2temp[,4]
  S2   <- s2temp[,3]
  s2ylim <- c(0,(1.1*max(s2temp[,4]+max(s2temp[,5]))))
  s2fit<- lm(s2res ~ S2)
  if (x == 1) {
    par(new=F)
    close.screen(1,2); split.screen(figs = c(1,2), erase=T)
    screen(1)
    plot(S1, s1res, las=1, cex=0.5, ylim=ylim, col=colo, ylab="Population mean petal number", xlab="Parameter 'select1'")
    abline(h=2.01, lwd=1.5, lty=2, col="blue")
    lines(xs1, predict(s1fit, list(S1 = seq(min(S1), max(S1), by=1))), col=colo)
    legend(130, magicnumber-((x-1)*(magicnumber/length(indices[,1]))), legend=paste("select2=",indices[x,1], sep=""),
           bty="n", pch=1, cex=0.65, pt.cex=0.5, col=colo)
    legend(50,0.1,  legend="Empirical pop. mean", pch=NA, lty=2, lwd=2, col="blue", bty="n", cex=0.65)
    screen(2)
    plot(S2, s2res, las=1, cex=0.5, ylim=ylim, col=colo, ylab="Population mean petal number", xlab="Parameter 'select2'")
    abline(h=2.01, lwd=1.5, lty=2, col="blue")
    lines(xs2, predict(s2fit, list(S2 = seq(min(S2), max(S2), by=0.01))), col=colo)
    legend(4.54, magicnumber-((x-1)*(magicnumber/length(indices[,1]))), legend=paste("select1=",indices[x,2], sep=""),
           bty="n", pch=1, cex=0.65, pt.cex=0.5, col=colo)
    legend(3.8,0.1,  legend="Empirical pop. mean", pch=NA, lty=2, lwd=2, col="blue", bty="n", cex=0.65)
  } else {
    screen(1, new=F)
    par(new=T)
    plot(S1, s1res, col=colo, axes=F, ylab="", xlab="", ylim=ylim, cex=0.5)
    lines(xs1, predict(s1fit, list(S1 = seq(min(S1), max(S1), by=1))), col=colo)
    legend(130, magicnumber-((x-1)*(magicnumber/length(indices[,1]))), legend=paste("select2=",indices[x,1], sep=""),
           bty="n", pch=1, cex=0.65, pt.cex=0.5, col=colo)
    screen(2, new=F)
    par(new=T)
    plot(S2, s2res, col=colo, ylim=ylim, axes=F, ylab="", xlab="", cex=0.5)
    lines(xs2, predict(s2fit, list(S2 = seq(min(S2), max(S2), by=0.01))), col=colo)
    legend(4.54, magicnumber-((x-1)*(magicnumber/length(indices[,1]))), legend=paste("select1=",indices[x,2], sep=""),
           bty="n", pch=1, cex=0.65, pt.cex=0.5, col=colo)
  }
}

dev.off()
