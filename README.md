# PredictionError

To install:
```R
devtools::install_github("matthewtyler/predictionError")
```

Simple Example:
```
n <- 2e3
n_v <- 150
n_t <- 100
n_p <- n - n_v - n_t
lab <- sample(n, n_v + n_t)
val <- sample(lab, n_v)
v <- 1 * (1:n %in% val)
p <- 1 * (! 1:n %in% lab)
t <- 1 - v - p

beta_true <- 2:4 / 10
sigma <- 5.0

Xu <- rnorm(n)
Xo <- rnorm(n)
epsilon <- sigma * rnorm(n)
y <- cbind(Xu, Xo, 1) %*% beta_true + epsilon
Zu <- Xu + 2.5 * rnorm(n) # Zu predicts Xu without being correlated with epsilon

Zu[t == 1] <- NA
Xu[p == 1] <- NA

missing_cov(y, Xu, Xo, Zu, v, t, p)
```