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
#' @param ER_test_signif_level default is 0.05; significance level for the ER test warning message, but note that the pvalue itself is not suppreseed
#' @param confint_signif_level default is 0.05; 1 - confint_signif_level determines the confidence level for the provided GMM estimator confidence intervals
#' @param ER_test default is TRUE; if TRUE, uses the validation sample to test the required exclusion restriction: `E(epsilon z_u) = 0' via Sargan's J-test results from the \code{gmm} package; see our paper for construction of the test
#' @param include_intercept default is TRUE; if TRUE, will append a columns of ones to Xo, inducing an intercept term in the model y ~ X
#' @param max_iter default is 25; determines how many iterations the GMM estimator will take before quitting
#' @param min_iter default is 2; min_iter - 1 determines how many iterations the GMM weighting matrix is based on both unlabeled and labeled data
#' @param tol default is 0.01; the algorithm uses [min_iter, max_iter] iterations, stopping if percent changes in beta are < tol for each element of beta
#' @param verbose default is TRUE; tells the function whether to print convergence and other warnings
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
  ER_test_signif_level = 0.05,
  confint_signif_level = 0.05,
  ER_test = TRUE, include_intercept = TRUE,
  min_iter = 2, max_iter = 25, tol = 0.01,
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

  Gamma_reg <- lm(X ~ Z - 1, subset = v == 1)
  # lm()'s Gamma needs to be transposed to be conformable with ours
  Gamma_hat <- t(coef(Gamma_reg))
  X_hat <- predict(Gamma_reg, newdata = data.frame(Z))
  eta_hat <- array(NA, dim = c(n, d_x))
  eta_hat[v == 1, ] <- resid(Gamma_reg)  

  beta_hat <- coef(lab_only <- lm(y ~ X - 1, subset = p == 0))

  B <- c(t(X[p == 0, ]) %*% y[p == 0] / (n_t + n_v),
    t(Z[t == 0, ]) %*% y[t == 0] / (n_p + n_v))
  A <- rbind(crossprod_self(X[p == 0, ]) / (n_t + n_v),
    t(Z[t == 0, ]) %*% X_hat[t == 0, ] / (n_p + n_v))

  lam_p <- n_v / (n_p + n_v)
  lam_t <- n_v / (n_t + n_v)

  new_beta <- function(W) {
    temp <- t(A) %*% W
    return(solve(temp %*% A) %*% temp %*% B)
  }

  Gamma_error <- array(0, dim = c(n, d_z, d_x))
  for (i in 1:n) {
    if (v[i] == 1) Gamma_error[i, , ] <- Z[i, ] %*% t(eta_hat[i, ])
  }

  weight_matrix <- function(beta) {

    m2_error <- array(0, dim = c(n, d_z))
    for (i in 1:n) {
      if (v[i] == 1) m2_error[i,] <- Gamma_error[i, , ] %*% beta
    }

    m1 <- array(0, dim = c(n, d_x))
    m2 <- array(0, dim = c(n, d_z))
    for (i in 1:n) {
      if (p[i] == 0) m1[i, ] <- lam_t * c(y[i] - X[i,] %*% beta) * X[i,]
      if (t[i] == 0) m2[i, ] <- lam_p * c(y[i] - X_hat[i, ] %*% beta) * Z[i,]
    }
    m2 <- m2 + m2_error
    m <- cbind(m1, m2)
    return(solve(crossprod_self(m) / n_v))
  }

  for (iter in 1:max_iter) {
    W <- weight_matrix(beta_hat)
    beta_old <- beta_hat
    beta_hat <- new_beta(W)
    if (iter >= min_iter & check_conv(beta_old, beta_hat, tol)) {
      if (verbose) message(paste0("Convergence after ", iter, " iterations!"))
      break
    }
  } 
  if (iter == max_iter & verbose) message(paste0("WARNING: Failure to converge after ", max_iter, "iterations!"))
  
  W <- weight_matrix(beta_hat)
  vcov <- solve(t(A) %*% W %*% A) / n_v
  se <- sqrt(diag(vcov))
  confint <- matrix(rep(beta_hat, 2), ncol = 2) +
    se %o% c(-1, 1) * qnorm(1 - confint_signif_level/2)
  
  if (ER_test) {
    pval <- perform_test(y[v == 1], X[v == 1,], all_instruments[v == 1,])
    if (verbose & pval < ER_test_signif_level) warning(paste0("WARNING: A statistical test of the required assumption `Exclusion Restriction: E(epsilon z_u) = 0' was rejected at the ", signif(ER_test_signif_level, 3), " significance level!"), immediate. = TRUE)
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

check_conv <- function(a, b, tol) {
  K <- length(a)
  for (k in 1:K)
    if (log(abs(a[k] - b[k])) - log(abs(a[k])) >= log(tol)) return(FALSE)
  return(TRUE)
}