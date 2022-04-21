# cohort size
#n <- c(100, 500, 1000, 2000)
n <- 100
cohorts <- 1000
data_generation <- function(cohorts = 1000, n, SEED = 1){
  set.seed(seed = SEED)
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
