# Bjorn Pieper, MPIPZ Cologne, Germany, July 2014
# simulated evolution of petal number in Cardamine hirsuta with and without selection
# mutations lead to a growing genepool by maintaining also the original lineage
# repeat the entire simulation a number of times
# probability to go extinct by selection against high petal number
# Here, selection is applied before truncating petal number to the range(0, 4)
# QTL effects are applied as if, their standard errors ignored int his version
# large number of simulations 
setwd("/home/pieper/MPIPZ/cluster/PN_evolve/run2_/production_run2")
##setwd("/home/bpuser/server/home/cluster/PN_evolve/run2_/")

empidat <- read.csv("/home/pieper/MPIPZ/x-perimentz/QTL_analysis/Meta/figures/acchisto/accpntn.csv", header = T)
##empidat <- read.csv("/home/bpuser/server/home/x-perimentz/QTL_analysis/Meta/figures/acchisto/accpntn.csv", header = T)

outfile <- "run2_v1_S2_4.1p_(new)_"
nrep <- 20000
save_every <- 9 
iter <- 1
apply_selection <- 0  # 0 is without selection (extinction), 1 is with. Switches after every save 

results <- matrix(nrow = nrep, ncol = 10)
resultz <- matrix(nrow = nrep, ncol = 45)
stats <- matrix(nrow = nrep, ncol = 6)  # colnames= c('MWU.twos','MWU.less','MWU.grea', 'KS.twos', 'KS.less','KS.grea'))
max_PN <- vector(mode = "numeric", length = nrep)
mu_PN <- vector(mode = "numeric", length = nrep)
min_PN <- vector(mode = "numeric", length = nrep)
PN4 <- vector(mode = "numeric", length = nrep)
PN0 <- vector(mode = "numeric", length = nrep)
cS1 <- vector(mode = "numeric", length = nrep)

select_no1 <- 148.04  # central value of S1 

calcP_ext <- function(x){
  return( 1/(1 + (select1^-(-select2 + (ancPN + sum(qtlevolve[, x]))))) )
}

init <- function() {
  nlineages <<- 1000  # number of lineages to model
  ancPN <<- 3.318  # ancestral petal number. When the sum of qtl effects is added a value of 2 petals is obtained 
  qtleffects <<- list(-0.386, c(-0.748, -0.438), 0.474, 0.626, c(-0.532, 0.490), c(0.538, 0.368), -0.442, -0.514, c(0.554, 0.464), 
                     -0.436, -0.912, -0.444, 0.534, c(0.922, 0.654), -0.544)
  qtlse <<- list(0.054, c(0.050, 0.065), 0.061, 0.059, c(0.057, 0.070), c(0.052, 0.049), 0.050, 0.051, c(0.050, 0.070),
                0.054, 0.063, 0.061, 0.065, c(0.048, 0.053), 0.067)
  p_mutate <<- 7.5e-04  # probability to mutate
#  select1 <<- 51.85  # this number determines the rate of increase in selection strength with increasing petal numbers
  select2 <<- 4.1  # this number determines at which petal number there is a 50% chance to go extinct 
  line_id <<- 1:nlineages
  linelife <<- rep(0, length(line_id))
  qtlevolve <<- matrix(nrow = 15, ncol = length(line_id))
  qtlevolve[, 1] <<- sapply(qtleffects, function(x){ sample(x,1) }) 
  linelife[1] <<- 1
  p_ext <<- rep(apply_selection, length(line_id))
  if (apply_selection == 1) { p_ext[1] <<- calcP_ext(1) } # calculate the probability to go extinct 
  panic <<- 0
}

exteval <- function() {
  return( ifelse(sum(linelife) > 1, runif(1, 0, 1), 1) )
}

drift <- function(j, k){
  if (length(qtleffects[[j]]) == 1) {
    return( ifelse(qtlevolve[j, new_life] != 0, 0, qtleffects[[j]]) )
  } else {
    return( sapply(list(c(0,qtleffects[[j]])), function(x){ sample(subset(x, x != qtlevolve[j, k]), 1) }) )
  }
}

extinction <- function(k) {
  qtlevolve[, k] <<- 0
  linelife[k] <<- 0
  p_ext[k] <<- 1
}
            
is.integer <- function(N) {
  !length(grep("[^[:digit:]]", as.character(N)))
}

  for (number in 1:nrep) {
    for ( vars1 in seq(0.90, 1.3, by=0.05) ){ 
    select1 <- select_no1 * vars1  
    init()
    while (sum(linelife) < nlineages) {
      for (k in subset(line_id, linelife == 1)) {
        if ( exteval() > p_ext[k]) {
          for (j in sample(1:length(qtleffects), length(qtleffects))) {
            if (runif(1, 0, 1) <= p_mutate) {
              if (sum(linelife) == nlineages) { break }
              new_life <- subset(line_id, linelife == 0)[1]
              qtlevolve[, new_life] <- qtlevolve[, k]
              linelife[new_life] <- 1
              qtlevolve[j, new_life] <- drift(j, k) 
              if (apply_selection == 1) {
                p_ext[new_life] <- calcP_ext(new_life) 
              }
            }
          }
        } else {
          extinction(k)
        }
      }
      panic = panic + 1
      if (panic == 50000) {
        print("PANIC!")
        break
      }
    }
    
    popPN <- vector(mode = "numeric", length = nlineages)
    
    for (x in 1:nlineages) {
      popPN[x] = ancPN + sum(qtlevolve[, x])
      if (popPN[x] > 4) {
        popPN[x] <- 4
      } else if (popPN[x] < 0) {
        popPN[x] <- 0
      }
    }
   
    lin_sample <- sample(1:nlineages, 45)
    resultz[iter, ] <- popPN[lin_sample]
    temp <- hist(popPN[lin_sample], breaks = seq(0, 4, by = 0.4), plot = F)
    results[iter, ] <- temp$counts
    max_PN[iter] <- max(popPN[lin_sample])
    mu_PN[iter] <- mean(popPN[lin_sample])
    min_PN[iter] <- min(popPN[lin_sample])
    PN4[iter] <- temp$counts[10]
    PN0[iter] <- temp$counts[1]
    cS1[iter] <- select1
    stats[[iter, 1]] <- wilcox.test(empidat[, 2], resultz[iter, ], alternative = "two.sided", 
      paired = FALSE, exact = F, correct = T, conf.int = F)[[3]]
    stats[[iter, 2]] <- wilcox.test(empidat[, 2], resultz[iter, ], alternative = "less", 
      paired = FALSE, exact = F, correct = T, conf.int = F)[[3]]
    stats[[iter, 3]] <- wilcox.test(empidat[, 2], resultz[iter, ], alternative = "greater", 
      paired = FALSE, exact = F, correct = T, conf.int = F)[[3]]
    stats[[iter, 4]] <- ks.test(empidat[, 2], resultz[iter, ], alternative = "two.sided", 
      paired = FALSE, exact = F, correct = T, conf.int = F)[[2]]
    stats[[iter, 5]] <- ks.test(empidat[, 2], resultz[iter, ], alternative = "less", 
      paired = FALSE, exact = F, correct = T, conf.int = F)[[2]]
    stats[[iter, 6]] <- ks.test(empidat[, 2], resultz[iter, ], alternative = "greater", 
      paired = FALSE, exact = F, correct = T, conf.int = F)[[2]]
    
    if (is.integer(iter/save_every)) {
      aux1 <- ifelse(iter == save_every, 1, iter - (save_every - 1))
      write.table(cbind(rep(apply_selection, save_every), cS1[aux1:iter], rep(select2, save_every), 
        aux1:iter, min_PN[aux1:iter], mu_PN[aux1:iter], max_PN[aux1:iter], 
        PN4[aux1:iter], PN0[aux1:iter]), file = paste(outfile, ".summ.tbl", sep=""), append = TRUE, 
        col.names = FALSE, row.names = FALSE)
      write.table(cbind(rep(apply_selection, save_every), cS1[aux1:iter], rep(select2, save_every), 
       stats[aux1:iter,]), file = paste(outfile, ".stats.tbl", sep=""), append = TRUE, 
        col.names = FALSE, row.names = FALSE)
      write.table(cbind(rep(apply_selection, save_every), cS1[aux1:iter], rep(select2, save_every), 
        results[aux1:iter,], resultz[aux1:iter,]), file = paste(outfile, ".resultsz.tbl", sep=""), append = TRUE, 
        col.names = FALSE, row.names = FALSE)
      write.table(cbind(rep(apply_selection, 15), 
        qtlevolve), file = paste(outfile, ".qtlevo.tbl", sep=""), append = TRUE, 
        col.names = FALSE, row.names = FALSE)
      apply_selection <- ifelse(apply_selection == 0, 1, 0) 
    } 
    iter <- iter + 1
  }
}
