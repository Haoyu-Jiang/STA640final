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
  for (i in 1:cohorts){
    res <- DR.estimator(X = as.matrix(mydata[[i]][, "X"]),
                            Y = as.matrix(mydata[[i]][, "Y"]),
                            Z = as.matrix(mydata[[i]][, c("Z1", "Z2", "Z3")]),
                            PS.Zid = PS.Zid,
                            OR.Zid = OR.Zid)
    taus[i] <- res[1]
    SE_ACMs[i] <- res[2]
  }
  Bias <- mean(taus)
  SE_ACM <- mean(SE_ACMs)
  SD <- sd(taus)
  
  return (c(Bias, SE_ACM, SD))
}

simu(cohorts = 1000, n = 1000, PS.Zid = c(1:3), OR.Zid = c(1, 3))

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

