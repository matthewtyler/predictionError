# predictionError

To install:
```R
devtools::install_github("matthewtyler/predictionError")
```

Simple Example:
```R
set.seed(1234)
n <- 2e3
n_v <- 150
n_t <- 100
n_p <- n - n_v - n_t
lab <- sample(n, n_v + n_t)
val <- sample(lab, n_v)
v <- 1 * (1:n %in% val)
p <- 1 * (! 1:n %in% lab)
t <- 1 - v - p

beta_true <- c(0.2, 0.4, 0.3)
sigma <- 1.0

Xu <- rnorm(n)
Xo <- rnorm(n)
epsilon <- sigma * rnorm(n)
y <- cbind(Xu, Xo, 1) %*% beta_true + epsilon
Zu <- Xu + rnorm(n) # Zu predicts Xu without being correlated with epsilon

Zu[t == 1] <- NA
Xu[p == 1] <- NA

predicted_covariates(y, Xu, Xo, Zu, v, t, p)
```