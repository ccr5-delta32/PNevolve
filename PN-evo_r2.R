# mutations lead to a growing genepool by maintaining also the original lineage
# repeat the entire simulation a number of times
# probability to go extinct by selection against high petal number
# The method is the same as above but no figures are plotted. Instead occurences of 4 petals,
# 0 petals, and max, mean and min petal number are scored with and without selection for a 
# large number of simulations 
setwd("/home/pieper/MPIPZ/cluster/PN_evolve")

is.integer <- function(N) {
  !length(grep("[^[:digit:]]", as.character(N)))
}

empi <- c(3, 1, 4, 5, 7, 8, 9, 6, 2, 0)
empidat <- read.csv("/home/pieper/MPIPZ/x-perimentz/QTL_analysis/Meta/figures/acchisto/accpntn.csv", 
  header = T)

para <- 1
nrep <- 2000
save_every <-10 
iter <- 1
apply_selection <- 0  # 0 is without selection (extinction), 1 is with 

results <- matrix(nrow = nrep, ncol = 10)
resultz <- matrix(nrow = nrep, ncol = 45)
stats <- matrix(nrow = nrep, ncol = 6)  # colnames= c('MWU.twos','MWU.less','MWU.grea', 'KS.twos', 'KS.less','KS.grea'))
max_PN <- vector(mode = "numeric", length = nrep)
mu_PN <- vector(mode = "numeric", length = nrep)
min_PN <- vector(mode = "numeric", length = nrep)
PN4 <- vector(mode = "numeric", length = nrep)
PN0 <- vector(mode = "numeric", length = nrep)

for (number in 1:nrep) {
  print(paste("repetition ", number, " out of ", nrep/2, " for selection = ", 
    apply_selection))
  nlineages <- 1000  # number of lineages to model
  ancPN <- 3.72  # ancestral petal number. When the sum of qtl effects is added a value of 2 petals is obtained 
  qtleffects <- c(-0.39, -0.75, 0.47, 0.63, -0.53, 0.54, -0.44, -0.51, 0.55, 
    -0.44, -0.91, -0.44, 0.53, 0.92, -0.54)
  p_mutate <- 7.5e-04  # probability to mutate
  select1 <- 75  # this number determines the rate of increase in selection strength with increasing petal numbers
  select2 <- 4.4  # this number determines at which petal number there is a 50% chance to go extinct 
  line_id <- 1:nlineages
  linelife <- rep(0, length(line_id))
  qtlevolve <- matrix(nrow = 15, ncol = length(line_id))
  qtlevolve[, 1] <- qtleffects
  linelife[1] <- 1
  if (apply_selection == 1) {
    p_ext <- rep(1, length(line_id))
    p_ext[1] <- 1/(1 + (select1^-(-select2 + (ancPN + sum(qtlevolve[, 1])))))  # calculate the probability to go extinct
  } else {
    p_ext <- rep(0, length(line_id))
  }
  panic <- 0
  
  while (sum(linelife) < nlineages) {
    for (k in subset(line_id, linelife == 1)) {
      if (sum(linelife) > 1) {
        temp1 <- runif(1, 0, 1)
      } else {
        temp1 <- 1
      }
      if (temp1 > p_ext[k]) {
        for (j in 1:length(qtleffects)) {
          if (runif(1, 0, 1) <= p_mutate) {
          if (sum(linelife) == nlineages) {
            break
          }
          new_life <- subset(line_id, linelife == 0)[1]
          qtlevolve[, new_life] <- qtlevolve[, k]
          linelife[new_life] <- 1
          if (qtlevolve[j, new_life] != 0) {
            qtlevolve[j, new_life] = 0
          } else {
            qtlevolve[j, new_life] = qtlevolve[j, k]
          }
          if (apply_selection == 1) {
            p_ext[new_life] <- 1/(1 + (select1^-(-select2 + (ancPN + 
            sum(qtlevolve[, new_life])))))  # calculate the probability to go extinct
          }
          break
          }
        }
      } else {
        qtlevolve[, k] <- 0
        linelife[k] <- 0
        p_ext[k] <- 1
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
  resultz[number, ] <- popPN[lin_sample]
  temp <- hist(popPN[lin_sample], breaks = seq(0, 4, by = 0.4), plot = F)
  results[number, ] <- temp$counts
  max_PN[number] <- max(popPN[lin_sample])
  mu_PN[number] <- mean(popPN[lin_sample])
  min_PN[number] <- min(popPN[lin_sample])
  PN4[number] <- temp$counts[10]
  PN0[number] <- temp$counts[1]
  stats[[number, 1]] <- wilcox.test(empidat[, 2], resultz[number, ], alternative = "two.sided", 
    paired = FALSE, exact = F, correct = T, conf.int = F)[[3]]
  stats[[number, 2]] <- wilcox.test(empidat[, 2], resultz[number, ], alternative = "less", 
    paired = FALSE, exact = F, correct = T, conf.int = F)[[3]]
  stats[[number, 3]] <- wilcox.test(empidat[, 2], resultz[number, ], alternative = "greater", 
    paired = FALSE, exact = F, correct = T, conf.int = F)[[3]]
  stats[[number, 4]] <- ks.test(empidat[, 2], resultz[number, ], alternative = "two.sided", 
    paired = FALSE, exact = F, correct = T, conf.int = F)[[2]]
  stats[[number, 5]] <- ks.test(empidat[, 2], resultz[number, ], alternative = "less", 
    paired = FALSE, exact = F, correct = T, conf.int = F)[[2]]
  stats[[number, 6]] <- ks.test(empidat[, 2], resultz[number, ], alternative = "greater", 
    paired = FALSE, exact = F, correct = T, conf.int = F)[[2]]
  
  if (is.integer(iter/save_every)) {
    if (iter == save_every) {
      aux1 = 1
    } else {
      aux1 = iter - (save_every - 1)
    }
    write.table(cbind(rep(para, save_every), rep(apply_selection, save_every), 
      aux1:iter, min_PN[aux1:iter], mu_PN[aux1:iter], max_PN[aux1:iter], 
      PN4[aux1:iter], PN0[aux1:iter]), file = "PN-evo_MWU+KS_200_4.5.summ.tbl", append = TRUE, 
      col.names = FALSE, row.names = FALSE)
    write.table(cbind(rep(para, save_every), rep(apply_selection, save_every), 
     stats[aux1:iter,]), file = "PN-evo_MWU+KS_200_4.5.stats.tbl", append = TRUE, 
      col.names = FALSE, row.names = FALSE)
    write.table(cbind(rep(para, save_every), rep(apply_selection, save_every), 
      results[aux1:iter,], resultz[aux1:iter,]), file = "PN-evo_MWU+KS_200_4.5.resultsz.tbl", append = TRUE, 
      col.names = FALSE, row.names = FALSE)
    write.table(cbind(rep(para, 15), rep(apply_selection, 15), 
      qtlevolve), file = "PN-evo_MWU+KS_200_4.5.qtlevo.tbl", append = TRUE, 
      col.names = FALSE, row.names = FALSE)
    if (apply_selection == 0) {
      apply_selection <- 1
    } else {
      apply_selection <- 0
    }
  }
  iter <- iter + 1
} 

