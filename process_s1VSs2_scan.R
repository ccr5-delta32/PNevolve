# Bjorn Pieper. MPIPZ, parCologne, Germany. Make pretty graphs of s1VSs2 scan results

setwd("/home/bpuser/server/home/cluster/PN_evolve/run2_")
scanmu1 <- read.table("run2_v1_s2VSs1_bis_.tbl", header=F)

## !!! uglyness alert !! ##
sc_V1 <- vector(mode="numeric", length=256)
sc_V2 <- vector(mode="numeric", length=256)
sc_V3 <- vector(mode="numeric", length=256)
sc_V4 <- vector(mode="numeric", length=256)
sc_V5 <- vector(mode="numeric", length=256)

scanmu <- cbind(sc_V1,sc_V2,sc_V3,sc_V4,sc_V5)
facepalm <- 1

for (x in unique(scanmu1[,2])){
  for (y in unique(scanmu1[,3])){
    scanmu[facepalm,1] <- facepalm
    scanmu[facepalm,2] <- mean(subset(scanmu1[,2], scanmu1[,2] == x))
    scanmu[facepalm,3] <- mean(subset(scanmu1[,3], scanmu1[,3] == y))
    scanmu[facepalm,4] <- mean(subset(scanmu1[,4], scanmu1[,2] == x & scanmu1[,3] == y))
    scanmu[facepalm,5] <- sd(subset(scanmu1[,4], scanmu1[,2] == x & scanmu1[,3] == y))/sqrt(length(subset(scanmu1[,4], scanmu1[,2] == x & scanmu1[,3] == y)))
    facepalm <- facepalm + 1
  }
}

xs1 <- seq(min(scanmu[,2]), max(scanmu[,2]), by=1)
xs2 <- seq(min(scanmu[,3]), max(scanmu[,3]), by=0.01)

indices <- cbind(seq(max(scanmu[,3]), min(scanmu[,3]), by=-0.1), seq(max(scanmu[,2]), min(scanmu[,2]), by=-10))

magicnumber <- 0.85
#png(file="scan_s1VSs2_results.png", width=3500, height=2000, unit="px", res=300)

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
  s1fit<- nls(s1res ~ (parA / S1) + (parB + log(S1, parC)), start=list(parA=-5, parB=1.3, parC=100))
  s2temp <- subset(scanmu, scanmu[,2] == as.character(indices[x,2]))
  s2res  <- s2temp[,4]
  S2   <- s2temp[,3]
  s2ylim <- c(0,(1.1*max(s2temp[,4]+max(s2temp[,5]))))
  s2fit<- lm(s2res ~ S2)
  print(indices[x,1])
  print(coef(s1fit))
  print(indices[x,2])
  print(coef(s2fit))
  if (x == 1) {
    par(new=F)
    close.screen(1,2); split.screen(figs = c(1,2), erase=T)
    screen(1)
    plot(S1, s1res, las=1, cex=0.5, ylim=ylim, col=colo, ylab="Population mean petal number", xlab="Parameter 'select1'")
    abline(h=2.01, lwd=1.5, lty=2, col="blue")
    lines(xs1, predict(s1fit, list(S1 = seq(min(S1), max(S1), by=1))), col=colo)
    legend(130, magicnumber-((x-1)*(magicnumber/length(indices[,1]))), legend=paste("select2=",indices[x,1], sep=""),
           bty="n", pch=1, cex=0.65, pt.cex=0.5, col=colo)
    legend(60,0.1,  legend="Empirical pop. mean", pch=NA, lty=2, lwd=2, col="blue", bty="n", cex=0.65)
    screen(2)
    plot(S2, s2res, las=1, cex=0.5, ylim=ylim, col=colo, ylab="Population mean petal number", xlab="Parameter 'select2'")
    abline(h=2.01, lwd=1.5, lty=2, col="blue")
    lines(xs2, predict(s2fit, list(S2 = seq(min(S2), max(S2), by=0.01))), col=colo)
    legend(4.54, magicnumber-((x-1)*(magicnumber/length(indices[,1]))), legend=paste("select1=",indices[x,2], sep=""),
           bty="n", pch=1, cex=0.65, pt.cex=0.5, col=colo)
    legend(3.9,0.1,  legend="Empirical pop. mean", pch=NA, lty=2, lwd=2, col="blue", bty="n", cex=0.65)
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

#dev.off()

# This function finds the proper value of S1 that leads to a popmu of 'des_mu' for a given S2
# take par_a, par_b, par_c from the output of the above loop
optiS1 <- function (par_a, par_b, par_c, des_mu) {
  P1 <- 15
  P2 <- 150
  accuracy <- 0.0000000001
    
  P3 <- P1 + ((P2 - P1) / 2)
  p3 <- (par_a / P3) + (par_b + (log(P3, par_c)))
  difference = abs(des_mu - p3)
  panic <- 0
  
  while(difference > accuracy) {
    panic <- panic + 1
    if (panic > 5000) { 
      print(paste("No solution found after", panic,"iterations"))
      break
    }
    
    p3 <- (par_a / P3) + (par_b + (log(P3, par_c)))

    difference = abs(des_mu - p3)

    if (p3 > des_mu) {
      P2 <- P3
      P3 <- P3 - ((P3 - P1) / 2)
    } else if (p3 < des_mu) {
      P1 <- P3
      P3 <- P3 + ((P2 - P3) / 2)
    }
  }
  print(paste("Solution found after",panic,"iterations:"))
  sprintf("%.16f",P3)
}

optiS1(-5.223803,1.385156,1939.959842,2.01)  #S2=4.1 / S1=148.0381
optiS1(-5.000568,1.392303,1023.960948,2.01)  #S2=4.2 / S1=101.7239
optiS1(-5.82006,1.67297,27581.78859 ,2.01)  #S2=4.3 / S1=71.8409
optiS1(-5.418639,1.646364,5050.985726,2.01)  #S2=4.4 / S1=53.0674
optiS1(-5.543990,1.775331,22366.393320,2.01)  #S2=4.5 / S1=40.8433
optiS1(-5.691476,1.858643,30332.379944,2.01)  #S2=4.6 / S1=31.2449
optiS1(-5.655510,1.893765,15796.529697,2.01)  #S2=4.7 / S1=25.7382
optiS1(-5.567027,2.028000,151665.267555,2.01)  #S2=4.8 / S1=20.5214

optiS1(-5.100110,1.454555,2315.848835,2.01)  #S2=4.1 / no solution
optiS1(-5.000568,1.392303,1023.960948,2.01)  #S2=4.2 / S1=106.98
optiS1(-5.381430,1.635917,3876.356114,2.01)  #S2=4.4 / S1=51.85
optiS1(-5.987016,1.956779,739289.039894,2.01)  #S2=4.6 / S1=30.12
optiS1(-5.813721,2.014233,84122.015062,2.01)  #S2=4.8 / S1=21.24"

optiS2 <- function(Int, slope, des_mu){
  result <- (des_mu - Int) / slope
  print(paste("The value of S2 that yields a population mu of", des_mu,"for the given S1 is:"))
  sprintf("%.16f",result)
}

mu_slope <- mean(c(0.6644838, 0.6310142, 0.66799, 0.6669765, 0.6631121, 0.6730596, 
  0.681241, 0.6890853, 0.6679688, 0.6718914, 0.6749862, 0.6921401, 0.7292495, 
  0.696658, 0.7078331, 0.5599476)) 

optiS2(-0.7504174,mu_slope,2.01)   #S1=155 / S2=4.1132585854176558
optiS2(-0.6066983,mu_slope,2.01)   #S1=145 / S2=3.8991048049192796
optiS2(-0.7665029,mu_slope,2.01)   #S1=135 / S2=4.1372273594790476
optiS2(-0.7710887,mu_slope,2.01)   #S1=125 / S2=4.1440605946343574
optiS2(-0.7735217,mu_slope,2.01)   #S1=115 / S2=4.1476859732232336
optiS2(-0.8261185,mu_slope,2.01)   #S1=105 / S2=4.2260597145152188
optiS2(-0.8772383,mu_slope,2.01)   #S1=95 / S2=4.3022325991792680
optiS2(-0.9259099,mu_slope,2.01)   #S1=85 / S2=4.3747574559512952
optiS2(-0.8672866,mu_slope,2.01)   #S1=75 / S2=4.2874037129881790
optiS2(-0.9171223,mu_slope,2.01)   #S1=65 / S2=4.3616631785622255
optiS2(-0.9627440,mu_slope,2.01)   #S1=55 / S2=4.4296434228565662
optiS2(-1.0748491,mu_slope,2.01)   #S1=45 / S2=4.5966896329855507
optiS2(-1.2889550,mu_slope,2.01)   #S1=35 / S2=4.9157257799695451
optiS2(-1.2670780,mu_slope,2.01)   #S1=25 / S2=4.4.8831271743843230
optiS2(-1.5240412,mu_slope,2.01)   #S1=15 / S2=5.2660243726617990
optiS2(-1.7038122,mu_slope,2.01)   #S1=5 / S2=5.5338985750049359

