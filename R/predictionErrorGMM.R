#' Generalized method of moments (GMM) estimator for overcoming
#' prediction error
#'
#' Combines ordinary least squares (OLS) on the hand-labeled set and
#' two-sample two-stage least squares (TS2SLS) on the unlabeled
#' set. See associated paper for details.
#'
#' @import stats
#' @import boot
#' @param Yl vector of outcome values for the mapping sample
#' @param Xl vector of X covariate values for the mapping sample
#' @param Zl matrix of Z predicted covariates values for the mapping sample
#' @param Vl matrix of V other covariate values for the mapping sample
#' @param Yu vector of outcome values for the unlabeled sample
#' @param Zu matrix of Z predicted covariates values for the unlabeled sample
#' @param Vu matrix of V other covariate values for the unlabeled sample
#' @param Yt vector of outcome values for the training sample (optional)
#' @param Xt matrix of X covariates values for the training sample (optional)
#' @param Vt matrix of V other covariate values for the training sample (optional)
#' @param num_boot number of bootstrap iterations to use in bootstrap asymptotic refinement; set to zero if you do not want to bootstrap the exclusion restriciton test p-value
#' @return A list of
#' \itemize{
#' \item{estimated coefficients \code{coefficients}}
#' \item{estimated standard errors \code{stderr}}
#' \item{covariance matrix estimate \code{vcov}}
#' \item{estimated coefficient of the residual on Z \code{rho_est} to evalute direction of exclusion restriction violation}
#' \item{exclusion restriction (overidentification) test statistic \code{OI_test_stat}}
#' \item{conventional chi-squared overidentificaiton test p-value \code{OI_test_pval_chisq}}
#' \item{bootstrap with asymptotic refinement overidentificaiton test p-value \code{OI_test_pval_boot}}
#' }
#' @export
predictionErrorGMM <- function(Yl, Xl, Zl, Vl = NULL,
                               Yu, Zu, Vu = NULL,
                               Yt = NULL, Xt = NULL, Vt = NULL, num_boot = 1e3) {
    mean_crossprod <- function(x, y) t(x) %*% y / nrow(x)

    to_boot_J_test <- function(XVl, indices, Yl, Zl, Rmean) {
        XVlb <- XVl[indices, ]
        Ylb <- Yl[indices]
        Zlb <- Zl[indices]
        theta_lb <- solve(t(XVlb) %*% XVlb) %*% t(XVlb) %*% Ylb
        rb <- sweep(c(Ylb - XVlb %*% theta_lb) * cbind(XVlb, Zlb), 2, Rmean)
        Rb <- colMeans(rb)
        Xib <- mean_crossprod(rb, rb)
        Qb <- nrow(XVlb) * t(Rb) %*% solve(Xib) %*% Rb
        return(c(Qb))
    }

    asymptotic_refinement <- function(XVl, Yl, Zl, R, Q, num_boot) {
        boot_res <- boot(XVl, to_boot_J_test, R = num_boot, Yl = Yl, Zl = Zl, Rmean = R)
        return(1 - mean(boot_res$t[, 1] < Q))
    }

    if (is.null(Vl)){
      Vl <- rep(1, length(Yl))
      Vu <- rep(1, length(Yu))
      if (!is.null(Yt)) Vt <- rep(1, length(Yt))
    }

    N_L <- length(Yl)
    N_U <- length(Yu)
    N_T <- length(Yt)
    N <- N_L + N_U + N_T
    lambda <- c(N_L, N_T, N_U)/N
    l <- c(rep(1, N_L), rep(0, N_T + N_U))
    t <- c(rep(0, N_L), rep(1, N_T), rep(0, N_U))
    u <- c(rep(0, N_L + N_T), rep(1, N_U))

    ZVl <- cbind(Zl, Vl)
    XVl <- cbind(Xl, Vl)
    ZVu <- cbind(Zu, Vu)
    ZVb <- rbind(ZVl, ZVu)
    XVt <- cbind(Xt, Vt)

    psi <- solve(t(ZVl) %*% ZVl) %*% t(ZVl) %*% Xl

    Xhatu <- ZVu %*% psi
    Xhatl <- ZVl %*% psi
    XhatVu <- cbind(Xhatu, Vu)
    XhatVl <- cbind(Xhatl, Vl)

    XVb <- rbind(XVl, XVt)
    theta <- solve(t(XVb) %*% XVb) %*% t(XVb) %*% c(Yl, Yt)
    theta_h <- theta
    beta <- theta[1]

    XhatVb <- rbind(XhatVl, XhatVu)

    B <- mean_crossprod(XhatVb, rbind(ZVl, ZVu)) %*% solve(mean_crossprod(ZVb, ZVb))

    g1 <- g2 <- array(0, c(N, 1 + ncol(Vl)))
    if (!is.null(Vt)) g1[t == 1, ] <- c(Yt - XVt %*% theta) * XVt
    g1[l == 1, ] <- c(Yl - XVl %*% theta) * XVl
    g2[l == 1 | u == 1, ] <- c(Yl - XhatVl %*% theta, Yu - XhatVu %*% theta) * XhatVb
    g2[l == 1, ] <- g2[l == 1, ] + (lambda[1] + lambda[3])/lambda[1] * beta *
        (as.vector(Xl - ZVl %*% psi) * ZVl) %*% t(B)

    g <- cbind(g1, g2)
    Omega <- mean_crossprod(g, g)
    W <- solve(Omega)
    G <-  -rbind(t(XVb) %*% XVb, t(XhatVb) %*% XhatVb) / N
    R1 <- rbind(t(XVb) %*% c(Yl, Yt), t(XhatVb) %*% c(Yl, Yu)) / N
    theta <- -solve(t(G) %*% W %*% G) %*% t(G) %*% W %*% R1
    beta <- theta[1]

    ### round 2

    g1 <- g2 <- array(0, dim = c(N, 1 + ncol(Vl)))
    if (!is.null(Vt)) g1[t == 1, ] <- c(Yt - XVt %*% theta) * XVt
    g1[l == 1, ] <- c(Yl - XVl %*% theta) * XVl
    g2[l == 1 | u == 1, ] <- c(Yl - XhatVl %*% theta, Yu - XhatVu %*% theta) * XhatVb
    g2[l == 1, ] <- g2[l == 1, ] + (lambda[1] + lambda[3])/lambda[1] * beta *
        (as.vector(Xl - ZVl %*% psi) * ZVl) %*% t(B)

    g <- cbind(g1, g2)
    Omega <- t(g) %*% g / N
    W <- solve(Omega)
    theta <- -solve(t(G) %*% W %*% G) %*% t(G) %*% W %*% R1
    beta <- theta[1]

    theta_var <- solve(t(G) %*% W %*% G) / N

    theta_l <- solve(t(XVl) %*% XVl) %*% t(XVl) %*% Yl
    r <- c(Yl - XVl %*% theta_l) * cbind(XVl, Zl)
    R <- colMeans(r)
    Xi <- mean_crossprod(r, r)
    Q <- c(N_L * t(R) %*% solve(Xi) %*% R)

    OI_test_stat <- Q
    OI_test_pval_chisq <- 1 - pchisq(Q, df = NCOL(Zl))
    OI_test_pval_boot <- asymptotic_refinement(XVl, Yl, Zl, R, Q, num_boot)

    beta_se <- sqrt(theta_var[1])

    eps_hat <- c(Yl - XVl %*% theta_h)
    rho_est <- coef(lm(eps_hat ~ Zl - 1))

    return(list(coef = theta,
                stderr = sqrt(diag(theta_var)),
                vcov = theta_var,
                rho_est = rho_est,
                OI_test_stat = OI_test_stat,
                OI_test_pval_chisq = OI_test_pval_chisq,
                OI_test_pval_boot = OI_test_pval_boot))
}
