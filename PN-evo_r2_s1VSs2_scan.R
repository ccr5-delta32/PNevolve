# Bjorn Pieper, MPIPZ Cologne, Germany, July 2014
# simulated evolution of petal number in Cardamine hirsuta with and without selection
# mutations lead to a growing genepool by maintaining also the original lineage
# repeat the entire simulation a number of times
# probability to go extinct by selection against high petal number
# Here, selection is applied before truncating petal number to the range(0, 4)
# QTL effects are applied as if, their standard errors ignored int his version
# THIS SCRIPT SCANS VALUES OF S1 AND S2 AND DETERMINES POPMU FOR EACH COMBINATION

## setwd("/home/pieper/MPIPZ/cluster/PN_evolve")
setwd("/home/bpuser/server/home/cluster/PN_evolve/run2_/")

# empidat <- read.csv("/home/pieper/MPIPZ/x-perimentz/QTL_analysis/Meta/figures/acchisto/accpntn.csv", header = T)
empidat <- read.csv("/home/bpuser/server/home/x-perimentz/QTL_analysis/Meta/figures/acchisto/accpntn.csv", header = T)

outfile <- "run2_v1_s2VSs1"
save_every <- 50 
nrep <- 16 * 16 * save_every 
iter <- 1
apply_selection <- 1  # 0 is without selection (extinction), 1 is with. Switches after every save 

mu_PN <- vector(mode = "numeric", length = nrep)
s1 <- vector(mode = "numeric", length = nrep)
s2 <- vector(mode = "numeric", length = nrep)
upopu <- vector(mode = "numeric", length = nrep)
s_mu <- vector(mode = "numeric", length = nrep/save_every)
s_se<- vector(mode = "numeric", length = nrep/save_every) 
s_s1 <- vector(mode = "numeric", length = nrep/save_every)
s_s2 <- vector(mode = "numeric", length = nrep/save_every) 

aux1 <- 1

init <- function() {
  nlineages <<- 1000  # number of lineages to model
  alpha <<- 3.318  # determined by 1000000 random combinations of 15 allelic effects
  qtleffects <<- list(-0.386, c(-0.748, -0.438), 0.474, 0.626, c(-0.532, 0.490), c(0.538, 0.368), -0.442, -0.514, c(0.554, 0.464), 
                     -0.436, -0.912, -0.444, 0.534, c(0.922, 0.654), -0.544)
  qtlse <<- list(0.054, c(0.050, 0.065), 0.061, 0.059, c(0.057, 0.070), c(0.052, 0.049), 0.050, 0.051, c(0.050, 0.070),
                0.054, 0.063, 0.061, 0.065, c(0.048, 0.053), 0.067)
  p_mutate <<- 7.5e-04  # probability to mutate
#  select1 <<- 75  # this number determines the rate of increase in selection strength with increasing petal numbers
#  select2 <<- 4.12  # this number determines at which petal number there is a 50% chance to go extinct 
  line_id <<- 1:nlineages
  linelife <<- rep(0, length(line_id))
  qtlevolve <<- matrix(nrow = 15, ncol = length(line_id))
  qtlevolve[, 1] <<- sapply(qtleffects, function(x){ sample(x,1) }) 
  linelife[1] <<- 1
  p_ext <<- rep(apply_selection, length(line_id))
  if (apply_selection == 1) { p_ext[1] <<- calcP_ext(1) } # calculate the probability to go extinct 
  panic <<- 0
}

calcP_ext <- function(x){
  return( 1/(1 + (select1^-(-select2 + (alpha + sum(qtlevolve[, x]))))) )
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

for (S2 in seq(3.3, 4.8, by=0.1)) {
  select2 <- S2
  for (S1 in seq(5, 155, by=10)) {
    select1 <- S1
    print(paste("select1:",S1," | select2:",S2))
    for (number in 1:save_every) {
      init()  
      while (sum(linelife) < nlineages) {
        for (k in subset(line_id, linelife == 1)) {
          if ( exteval() > p_ext[k]) {
            for (j in 1:length(qtleffects)) {
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
        popPN[x] = alpha + sum(qtlevolve[, x])
        if (popPN[x] > 4) {
          popPN[x] <- 4
        } else if (popPN[x] < 0) {
          popPN[x] <- 0
        }
      }
     
      lin_sample <- sample(1:nlineages, 45)
      mu_PN[iter] <- mean(popPN[lin_sample])
      s1[iter] <- select1
      s2[iter] <- select2
      iter = iter +1
    }
    s_mu[(iter-1)/save_every] <- mean(mu_PN[aux1:(iter-1)])
    s_se[(iter-1)/save_every] <- sd(mu_PN[aux1:(iter-1)])/sqrt(length(mu_PN[aux1:(iter-1)]))
    s_s1[(iter-1)/save_every] <- select1
    s_s2[(iter-1)/save_every] <- select2

    write.table(cbind(aux1:(iter-1), s1[aux1:(iter-1)], s2[aux1:(iter-1)], mu_PN[aux1:(iter-1)]),
      file = paste(outfile, ".tbl", sep=""), append = TRUE, col.names = FALSE,
      row.names = FALSE)
    aux1 = aux1 + save_every
  }
}

write.table(cbind(1:nrep/save_every, s_s1, s_s2, s_mu, s_se),
  file = paste(outfile, ".mu.tbl", sep=""), append = TRUE, col.names = FALSE,
  row.names = FALSE)

