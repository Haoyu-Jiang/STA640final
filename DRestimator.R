# Doubly Robust Estimator
# X: the binary (0-1) treatment of interest! dimension n*1
# Y: binary outcome (0-1), dimension n*1
# Z: confounders, dimension n*p, all of them are continuous or 0-1 binary
# For the potential outcome models, we refered to codes provided by Prof. Fan Li
DR.estimator <- function(X, Y, Z){
  # use logistic regression to estimate propensity score
  pscore <- glm(X ~ Z, family = "binomial")$fitted.values
  
  # use logistic regression to predict potential outcomes
  m1 <- glm(Y ~ Z, weights = X, family = "gaussian")$fitted.values
  m0 <- glm(Y ~ Z, weights = (1 - X), family = "gaussian")$fitted.values
  
  # doubly robost estimator
  tau.hat <- mean(m1 + X * (Y - m1) / pscore) - mean(m0 + (1 - X) * (Y - m0) / (1 - pscore))
  
  return (tau.hat)
}

# simple example
set.seed(123)
n <- 10000
Z <- matrix(rnorm(n * 3), nrow = n)
beta.ps <- as.matrix(c(1, 2, 3))
ps <- exp(Z %*% beta.ps) / (1 + exp(Z %*% beta.ps))
X <- rbinom(n, 1, prob = ps)

beta.Y <- as.matrix(c(0, 1, 0, -2)) # true causal effect = 0
Y.prob <- exp(cbind(X, Z) %*% beta.Y) / (1 + exp(cbind(X, Z) %*% beta.Y))
Y <- rbinom(n, 1, prob = Y.prob)

DR.estimator(X, Y, Z)



