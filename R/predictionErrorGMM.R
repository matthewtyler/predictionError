#' GMM estimator for models with (imperfectly) predicted covariates
#'
#' Optimally combines OLS and 2SLS on labeled and unlabeled data,
#' given an exclusion restriction. See associated paper for details.
#'
#' @import stats
#' @import gmm
#' @param y vector of n outcome values
#' @param Xu matrix/vector of n possibly unobserved (i.e., NA) covariate values; must be observed when v or t == 1
#' @param Xo matrix/vector of n fully observed covariate values; may be set to NULL, but should never contain a constant term --- control the inclusion of a constant term via \code{include_intercept} option
#' @param Zu matrix/vector of n fully observed predicted Xu values; must be observed when v or p == 1
#' @param v vector of n 1/0 where v[i] == 1 if unit i is validation sample, == 0 otherwise; sum(v) > 0 is required
#' @param t vector of n 1/0 where v[i] == 1 if unit i is training sample, == 0 otherwise; sum(t) == 0 is allowed
#' @param p vector of n 1/0 where v[i] == 1 if unit i is primary sample, == 0 otherwise; sum(p) == 0 is allowed, but then it wouldn't make sense to use this package
#' @param max_iter default is 1 for two-step GMM; the GMM estimator estimates beta max_iter + 1 times, with the first step using the labeled-only estimator to calculate the weighting matrix
#' @param ER_test_signif_level default is 0.05; significance level for the ER test warning message, but note that the pvalue itself is not suppreseed
#' @param confint_signif_level default is 0.05; 1 - confint_signif_level determines the confidence level for the provided GMM estimator confidence intervals
#' @param ER_test default is TRUE; if TRUE, uses the validation sample to test the required exclusion restriction: `E(epsilon z_u) = 0' via Sargan's J-test results from the \code{gmm} package; see our paper for construction of the test
#' @param include_intercept default is TRUE; if TRUE, will append a columns of ones to Xo, inducing an intercept term in the model y ~ X
#' @return A list of
#' \itemize{
#' \item{GMM estimated coefficients \code{beta}}
#' \item{estimated standard errors \code{se}}
#' \item{estimated confidence intervals \code{confint} at the confidence level 1 - \code{confint_signif_level}}
#' \item{variance-covariance matrix estimate \code{vcov}}
#' \item{exclusion restriction test p-value \code{pval}}
#' \item{labeled-only estimator of beta \code{lab_only}}
#' }
#' @examples
#' set.seed(1234)
#' n <- 2e3
#' n_v <- 150
#' n_t <- 100
#' n_p <- n - n_v - n_t
#' lab <- sample(n, n_v + n_t)
#' val <- sample(lab, n_v)
#' v <- 1 * (1:n %in% val)
#' p <- 1 * (! 1:n %in% lab)
#' t <- 1 - v - p
#' 
#' beta_true <- c(0.2, 0.4, 0.3)
#' sigma <- 1.0
#' 
#' Xu <- rnorm(n)
#' Xo <- rnorm(n)
#' epsilon <- sigma * rnorm(n)
#' y <- cbind(Xu, Xo, 1) %*% beta_true + epsilon
#' Zu <- Xu + rnorm(n) # Zu predicts Xu without being correlated with epsilon
#' 
#' Zu[t == 1] <- NA
#' Xu[p == 1] <- NA
#' 
#' predicted_covariates(y, Xu, Xo, Zu, v, t, p)
#' @export
predicted_covariates <- function(y, Xu, Xo, Zu, v, t, p,
  max_iter = 1, ER_test_signif_level = 0.05,
  confint_signif_level = 0.05,
  ER_test = TRUE, include_intercept = TRUE,
  verbose = TRUE) {
  
  n_v <- sum(v)
  n_t <- sum(t)
  n_p <- sum(p)
  n <- NROW(y)

  if (include_intercept) {
    if (verbose) message("Appended column of ones to Xo. If that is undesired, re-run function with include_intercept = FALSE.")
    if (is.null(Xo)) {
      Xo <- rep(1, n)
    } else {
      Xo <- cbind(Xo, 1) 
    }
  } else {
    if (is.null(Xo)) {
      stop("Sorry! Package currently does not support models without any Xo and without an intercept.")
    }
  }

  X <- cbind(Xu, Xo)
  Z <- cbind(Zu, Xo)
  all_instruments <- cbind(X, Zu)

  d_x <- NCOL(X)
  d_z <- NCOL(Z)
  d_xz <- d_x * d_z

  alpha_reg <- lm(X ~ Z - 1, subset = v == 1)
  alpha_hat <- matrix(coef(alpha_reg), ncol = 1)
  eta_hat <- array(NA, dim = c(n, d_x))
  eta_hat[v == 1, ] <- resid(alpha_reg)
  X_hat <- predict(alpha_reg, newdata = data.frame(Z))

  beta_hat <- coef(lab_only <- lm(y ~ X - 1, subset = p == 0))

  b <- c(t(X[p == 0, ]) %*% y[p == 0] / (n_t + n_v),
    t(X_hat[t == 0, ]) %*% y[t == 0] / (n_p + n_v))
  A <- rbind(crossprod_self(X[p == 0, ]) / (n_t + n_v),
    crossprod_self(X_hat[t == 0, ]) / (n_p + n_v))

  lam_p <- n_v / (n_p + n_v)
  lam_t <- n_v / (n_t + n_v)

  new_beta <- function(W) {
    temp <- t(A) %*% W
    return(solve(temp %*% A) %*% temp %*% b)
  }

  alpha_error <- array(0, dim = c(n, d_xz))
  for (i in 1:n) {
    if (v[i] == 1) {
      alpha_error[i, ] <- Z[i, ] * rep(eta_hat[i, ], each = d_x)
      ## above is equiavlent to t(I_{d_x} \otimes Z[i,]) %*% eta_hat[i, ]  
    }
  }

  L2 <- crossprod_self(Z[t == 0, ]) / (n_p + n_v)
  L2 <- kronecker(diag(d_x), L2)
  L2 <- solve(L2)

  weight_matrix <- function(beta) {

    L1 <- array(0, dim = c(n, d_x, d_xz))
    for (i in 1:n) {
      if (v[i] + p[i] > 0) L1[i, , ] <- X_hat[i, ] %*% t(rep(beta, each = d_x) * Z[i,])
      ## above is equivalent to t(beta) %*% (I_{d_x} \otimes Z[i,])
    }
    L1 <- apply(L1, c(2, 3), sum) / (n_p + n_v)

    L <- L1 %*% L2
    alpha_error <- alpha_error %*% t(L) # same as left-multiplying the rows by L

    m1 <- array(0, dim = c(n, d_x))
    m2 <- array(0, dim = c(n, d_x))
    for (i in 1:n) {
      if (p[i] == 0) m1[i, ] <- lam_t * c(y[i] - X[i,] %*% beta) * X[i,]
      if (t[i] == 0) m2[i, ] <- lam_p * c(y[i] - X_hat[i, ] %*% beta) * X_hat[i,]
    }
    m2 <- m2 + alpha_error
    m <- cbind(m1, m2)
    return(solve(crossprod_self(m) / n_v))
  }

  for (iter in 1:max_iter) {
    W <- weight_matrix(beta_hat)
    beta_hat <- new_beta(W)
  } 
  
  W <- weight_matrix(beta_hat)
  vcov <- solve(t(A) %*% W %*% A) / n_v
  se <- sqrt(diag(vcov))
  confint <- matrix(rep(beta_hat, 2), ncol = 2) +
    se %o% c(-1, 1) * qnorm(1 - confint_signif_level/2)
  
  if (ER_test) {
    pval <- perform_test(y[v == 1], X[v == 1,], all_instruments[v == 1,])
    if (verbose & pval < ER_test_signif_level) message(paste0("WARNING: A statistical test of the required assumption `Exclusion Restriction: E(epsilon z_u) = 0' was rejected at the ", signif(ER_test_signif_level, 3), " significance level!"))
  } else {
    pval <- NA
  }

  return(list(beta = beta_hat, se = se, confint = confint,
    vcov = vcov, ER_pval = pval, lab_only = coef(lab_only)))
}

crossprod_self <- function(X) t(X) %*% X

perform_test <- function(y, X, Z) {
  res <- gmm(g = y ~ X - 1, x = Z, vcov = "iid")
  pval <- summary(res)$stest$test[2]
  return(pval)
}