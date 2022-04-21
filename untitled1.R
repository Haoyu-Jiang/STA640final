library(tidyverse)

mydata <- data_generation(cohorts = 1000, n = 2000, SEED = 1)

taus <- rep(NA, length(mydata))
for (i in 1:length(mydata)){
  taus[i] <- DR.estimator(X = as.matrix(mydata[[i]][, "X"]),
                          Y = as.matrix(mydata[[i]][, "Y"]),
                          Z = as.matrix(mydata[[i]][, c("Z1", "Z2", "Z3")]),
                          PS.Zid = c(1:3),
                          OR.Zid = c(1))
}
mean(taus)
sd(taus)

mydata[[1]] %>% pull(X) %>% mean()

tmp <- rep(NA, 1000)
for(i in 1:1000){
  Y <- mydata[[i]] %>% pull(Y)
  X <- mydata[[i]] %>% pull(X)
  tmp[i] <- mean(Y[X == 1]) - mean(Y[X == 0])
}
mean(tmp)
















