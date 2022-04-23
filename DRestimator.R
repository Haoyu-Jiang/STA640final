data_generation <- function(cohorts = 1000, n){
  data <- list()
  for (i in 1:cohorts){
    Z1 <- rnorm(n, mean = 0, sd = 1)
    Z2 <- rnorm(n, mean = 0, sd = 1)
    Z3 <- rbinom(n, 1, prob = 0.3)
    # we use true model to generate x
    beta0 <- 1.5
    beta1 <- 1
    beta2 <- -2
    beta3 <- 1
    p <- 1 / (1 + exp(-(beta0 + beta1 * Z1 + beta2 * Z2 + beta3 * Z3)))
    
    # X <- as.integer(p + runif(n) < 0.91)
    X <- rbinom(n, 1, prob = p)
    # then generate Y
    epsilon <- rnorm(n, 0, sd = sqrt(0.08))
    Y <- Z1 + Z3 + 2 * rnorm(n, 0, 1) + epsilon
    data[[i]] <- data.frame(X = X,
                            Z1 = Z1,
                            Z2 = Z2,
                            Z3 = Z3,
                            Y = Y)
  }
  
  return (data)
}


# Doubly Robust Estimator
# X: the binary (0-1) treatment of interest! dimension n*1
# Y: continuous outcome, dimension n*1
# Z: confounders, dimension n*p, all of them are continuous or 0-1 binary
# PS.Zid: index of confounders used in PS model
# OR.Zid: index of confounders used in OR model.
# For the potential outcome models, we refered to codes provided by Prof. Fan Li
DR.estimator <- function(X, Y, Z, PS.Zid, OR.Zid){
  # use logistic regression to estimate propensity score
  pscore <- glm(X ~ Z[, PS.Zid], family = "binomial")$fitted.values
  
  # use logistic regression to predict potential outcomes
  m1 <- glm(Y ~ Z[, OR.Zid], weights = X, family = "gaussian")$fitted.values
  m0 <- glm(Y ~ Z[, OR.Zid], weights = (1 - X), family = "gaussian")$fitted.values
  
  # doubly robost estimator
  tau.hat <- mean(X * Y / pscore - m1 * (X - pscore) / pscore) -
    mean((1 - X) * Y / (1 - pscore) + m0 * (X - pscore) / (1 - pscore))
  
  tauis <- (X * Y / pscore - m1 * (X - pscore) / pscore) - 
    ((1 - X) * Y / (1 - pscore) + m0 * (X - pscore) / (1 - pscore))
  seACM <- sqrt(sum((tauis - tau.hat)^2) / length(Y)^2)

  return (c(tau.hat, seACM))
}


simu <- function(cohorts = 1000, n, PS.Zid, OR.Zid){
  mydata <- data_generation(cohorts = cohorts, n = n)
  
  # taus: store 1000 tau.dr estimates at each cohort
  taus <- rep(NA, cohorts)
  SE_ACMs <- rep(NA, cohorts)
  percentile.cover <- rep(NA, cohorts)
  SEACM.cover <- rep(NA, cohorts)
  SESD.cover <- rep(NA, cohorts)
  
  for (i in 1:cohorts){
    if (i %% 100 == 0){
      print(i)
    }
    res <- DR.estimator(X = as.matrix(mydata[[i]][, "X"]),
                            Y = as.matrix(mydata[[i]][, "Y"]),
                            Z = as.matrix(mydata[[i]][, c("Z1", "Z2", "Z3")]),
                            PS.Zid = PS.Zid,
                            OR.Zid = OR.Zid)
    taus[i] <- res[1]
    SE_ACMs[i] <- res[2]
    
    bsize <- 1000 # bootstrap times
    taus.boot <- rep(NA, bsize)
    for (j in 1:bsize){
      idx <- sample(n, replace = TRUE)
      res.boot <- DR.estimator(X = as.matrix(mydata[[i]][idx, "X"]),
                          Y = as.matrix(mydata[[i]][idx, "Y"]),
                          Z = as.matrix(mydata[[i]][idx, c("Z1", "Z2", "Z3")]),
                          PS.Zid = PS.Zid,
                          OR.Zid = OR.Zid)
      taus.boot[j] <- res.boot[1]
    }
    
    CI.percentile <- quantile(taus.boot, c(0.025, 0.975))
    percentile.cover[i] <- (CI.percentile[1] < 0 & 0 < CI.percentile[2])
    
    SEACM.cover[i] <- (taus[i] - 1.96 * SE_ACMs[i] < 0 & 0 < taus[i] + 1.96 * SE_ACMs[i])
    SESD.cover[i] <- (taus[i] - 1.96 * sd(taus.boot)< 0 & 0 < taus[i] + 1.96 * sd(taus.boot))
    
  }
  Bias <- mean(taus)
  SE_ACM <- mean(SE_ACMs)
  SD <- sd(taus)
  ratio = SE_ACM / SD
  SE_ACM_Coverage <- mean(SEACM.cover)
  SE_SD_Coverage <- mean(SESD.cover)
  Percentile_Coverage <- mean(percentile.cover)
  SE_ACM_Coverage_CI <- c((2 * cohorts * SE_ACM_Coverage + 1.96^2 - 1.96 * sqrt(1.96^2 + 4 * cohorts * SE_ACM_Coverage * (1 - SE_ACM_Coverage))) / (2 * (cohorts + 1.96^2)), 
                          (2 * cohorts * SE_ACM_Coverage + 1.96^2 + 1.96 * sqrt(1.96^2 + 4 * cohorts * SE_ACM_Coverage * (1 - SE_ACM_Coverage))) / (2 * (cohorts + 1.96^2)))
  SE_SD_Coverage_CI <- c((2 * cohorts * SE_SD_Coverage + 1.96^2 - 1.96 * sqrt(1.96^2 + 4 * cohorts * SE_SD_Coverage * (1 - SE_SD_Coverage))) / (2 * (cohorts + 1.96^2)), 
                          (2 * cohorts * SE_SD_Coverage + 1.96^2 + 1.96 * sqrt(1.96^2 + 4 * cohorts * SE_SD_Coverage * (1 - SE_SD_Coverage))) / (2 * (cohorts + 1.96^2)))
  Percentile_Coverage_CI <- c((2 * cohorts * Percentile_Coverage + 1.96^2 - 1.96 * sqrt(1.96^2 + 4 * cohorts * Percentile_Coverage * (1 - Percentile_Coverage))) / (2 * (cohorts + 1.96^2)), 
                          (2 * cohorts * Percentile_Coverage + 1.96^2 + 1.96 * sqrt(1.96^2 + 4 * cohorts * Percentile_Coverage * (1 - Percentile_Coverage))) / (2 * (cohorts + 1.96^2)))
  
  return (c(Bias, SE_ACM, SD, ratio, SE_ACM_Coverage, SE_ACM_Coverage_CI, SE_SD_Coverage, SE_SD_Coverage_CI, Percentile_Coverage, Percentile_Coverage_CI))
}

ttt <- simu(cohorts = 1000, n = 100, PS.Zid = c(1:3), OR.Zid = c(1, 3))

# scenario 1
res1.df <- data.frame()
for (n in c(100, 500, 1000, 2000)){
  res1.df <- rbind(
    res1.df,
    c(n, simu(cohorts = 1000, n = n, PS.Zid = c(1:3), OR.Zid = c(1, 3)))
  )
}

save(res1.df, file = "res1.df.RData")

# scenario 2
res2.df <- data.frame()
for (n in c(100, 500, 1000, 2000)){
  res2.df <- rbind(
    res2.df,
    c(n, simu(cohorts = 1000, n = n, PS.Zid = c(1), OR.Zid = c(1, 3)))
  )
}

# scenario 3
res3.df <- data.frame()
for (n in c(100, 500, 1000, 2000)){
  res3.df <- rbind(
    res3.df,
    c(n, simu(cohorts = 1000, n = n, PS.Zid = c(1:3), OR.Zid = c(1)))
  )
}
save(res3.df, file = "res3.df.RData")

# simple example
set.seed(12345)

tmp <- rep(NA, 1000)
for (i in 1:1000){
  n <- 1000
  Z <- matrix(rnorm(n * 3), nrow = n)
  beta.ps <- as.matrix(c(1, 2, 3))
  ps <- exp(Z %*% beta.ps) / (1 + exp(Z %*% beta.ps))
  X <- rbinom(n, 1, prob = ps)
  
  beta.Y <- as.matrix(c(3, 1, 0, -2)) # true causal effect = 3
  Y <- cbind(X, Z) %*% beta.Y + 0.5 * rnorm(n)
  
  tmp[i] <- DR.estimator(X, Y, Z, PS.Zid = c(1:3), OR.Zid = c(1, 3))[1]
}

mean(tmp)


df <- data.frame(Y = Y , Z1 = Z[, 1], Z2 = Z[, 2], Z3 = Z[, 3], X = X)
glm1 <- glm(Y ~ Z1 + Z2 + Z3, data = df %>% filter(X == 0))
m00 <- predict(glm1, newdata = df)
cbind(m00, m0)

