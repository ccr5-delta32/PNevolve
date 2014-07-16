qtleffects <<- list(-0.386, c(-0.748, -0.438), 0.474, 0.626, c(-0.532, 0.490), c(0.538, 0.368), -0.442, -0.514, c(0.554, 0.464), 
                    -0.436, -0.912, -0.444, 0.534, c(0.922, 0.654), -0.544)

nsim <- 1000000
alpha <- vector(mode="numeric", length=nsim)

for (sim in 1:nsim) {
  temp <- vector(mode="numeric", length=15)
  for (x in 1:15) {
    if (length(qtleffects[[x]]) == 1) {
      temp[x] <- qtleffects[[x]]
    } else {
      temp[x] <- sample(qtleffects[[x]],1)
    }
  alpha[sim] <- 2.41 - sum(temp)
  }
}

ualpha <- mean(alpha)
print(ualpha) 
