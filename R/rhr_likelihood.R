# Bivariate Left-Censored Log-Likelihood and Posterior Evaluation

#' Bivariate Left-Censored Log-Likelihood under Reversed Hazard Rate
#'
#' Evaluates the log-likelihood function for bivariate survival data with left censoring
#' and shared frailty based on reversed hazard rates.
#'
#' @param params Numeric vector of all model parameters (frailty parameters, baseline parameters for T1 and T2, and regression coefficients).
#' @param t1 Numeric vector of first failure/censoring times.
#' @param t2 Numeric vector of second failure/censoring times.
#' @param d1 Numeric or integer vector of censoring indicators for first component (1 = observed, 0 = left-censored).
#' @param d2 Numeric or integer vector of censoring indicators for second component (1 = observed, 0 = left-censored).
#' @param X Numeric matrix of covariates.
#' @param frailty Character string specifying frailty distribution.
#' @param baseline Character string specifying baseline distribution.
#' @param separate_baseline Logical indicating whether baseline parameters for T1 and T2 are estimated separately. Default is \code{TRUE}.
#'
#' @return Numeric scalar representing total log-likelihood.
#' @export
#'
#' @examples
#' data(austrian_twin)
#' params <- c(zeta1 = 12.4, delta1 = 0.11, zeta2 = 11.8, delta2 = 0.11,
#'             theta = 3.4, beta1 = 0.007, beta2 = -0.11, beta3 = 0.13)
#' X <- as.matrix(austrian_twin[, c("age", "k2", "k3")])
#' loglik_biv_rhr(params, austrian_twin$t1, austrian_twin$t2,
#'                austrian_twin$d1, austrian_twin$d2,
#'                X = X, frailty = "WL", baseline = "GE")
loglik_biv_rhr <- function(params, t1, t2, d1, d2, X = NULL, frailty = c("WL", "GL", "GIG", "Gamma", "IG", "Lindley", "None"), baseline = c("GE", "EG", "GIR", "MIW", "GLL-I", "GLL-II", "GIE", "Weibull", "Exponential", "Log-Logistic"), separate_baseline = TRUE) {
  frailty <- match.arg(frailty)
  baseline <- match.arg(baseline)
  
  b_names <- switch(baseline,
    "GE" = c("zeta", "delta"),
    "EG" = c("zeta", "delta"),
    "GIR" = c("zeta", "delta"),
    "MIW" = c("zeta", "delta", "xi"),
    "GLL-I" = c("lambda", "alpha", "gamma"),
    "GLL-II" = c("lambda", "alpha", "gamma"),
    "GIE" = c("zeta", "delta"),
    "Weibull" = c("zeta", "delta"),
    "Exponential" = c("delta"),
    "Log-Logistic" = c("lambda", "gamma")
  )
  n_b <- length(b_names)
  
  idx <- 1
  if (separate_baseline) {
    p1_vals <- params[idx:(idx + n_b - 1)]
    names(p1_vals) <- b_names
    idx <- idx + n_b
    p2_vals <- params[idx:(idx + n_b - 1)]
    names(p2_vals) <- b_names
    idx <- idx + n_b
  } else {
    p1_vals <- params[idx:(idx + n_b - 1)]
    names(p1_vals) <- b_names
    p2_vals <- p1_vals
    idx <- idx + n_b
  }
  
  if (frailty == "None") {
    theta <- 1.0
    mu <- NULL
    theta2 <- NULL
  } else if (frailty == "GL") {
    theta <- as.numeric(params[idx])
    idx <- idx + 1
    mu <- as.numeric(params[idx])
    idx <- idx + 1
    theta2 <- NULL
  } else if (frailty == "GIG") {
    theta <- as.numeric(params[idx])
    idx <- idx + 1
    theta2 <- as.numeric(params[idx])
    idx <- idx + 1
    mu <- NULL
  } else {
    theta <- as.numeric(params[idx])
    idx <- idx + 1
    mu <- NULL
    theta2 <- NULL
  }
  
  p_cov <- if (!is.null(X)) ncol(as.matrix(X)) else 0
  if (p_cov > 0 && idx <= length(params)) {
    beta <- params[idx:length(params)]
  } else {
    beta <- NULL
  }
  
  if (any(p1_vals <= 0) || any(p2_vals <= 0)) return(-1e10)
  if (frailty != "None" && theta <= 0) return(-1e10)
  if (frailty == "GL" && (mu <= 0 || mu >= (1 + theta))) return(-1e10)
  if (frailty == "GIG" && theta2 <= 0) return(-1e10)
  
  n <- length(t1)
  if (!is.null(X) && !is.null(beta)) {
    eta <- exp(as.numeric(as.matrix(X) %*% beta))
    eta1 <- eta
    eta2 <- eta
  } else {
    eta1 <- rep(1, n)
    eta2 <- rep(1, n)
  }
  
  r01 <- rhr_baseline(t1, baseline, p1_vals)
  r02 <- rhr_baseline(t2, baseline, p2_vals)
  R01 <- crhr_baseline(t1, baseline, p1_vals)
  R02 <- crhr_baseline(t2, baseline, p2_vals)
  
  s <- R01 * eta1 + R02 * eta2
  
  L <- laplace_frailty(s, frailty, theta, mu, theta2)
  dL <- dlaplace_frailty(s, frailty, theta, mu, theta2)
  d2L <- d2laplace_frailty(s, frailty, theta, mu, theta2)
  
  log_L <- numeric(n)
  
  c1 <- which(d1 == 1 & d2 == 1)
  if (length(c1) > 0) {
    log_L[c1] <- log(r01[c1]) + log(eta1[c1]) + log(r02[c1]) + log(eta2[c1]) + log(pmax(d2L[c1], 1e-15))
  }
  
  c2 <- which(d1 == 1 & d2 == 0)
  if (length(c2) > 0) {
    log_L[c2] <- log(r01[c2]) + log(eta1[c2]) + log(pmax(-dL[c2], 1e-15))
  }
  
  c3 <- which(d1 == 0 & d2 == 1)
  if (length(c3) > 0) {
    log_L[c3] <- log(r02[c3]) + log(eta2[c3]) + log(pmax(-dL[c3], 1e-15))
  }
  
  c4 <- which(d1 == 0 & d2 == 0)
  if (length(c4) > 0) {
    log_L[c4] <- log(pmax(L[c4], 1e-15))
  }
  
  log_L[!is.finite(log_L)] <- -1e10
  
  return(sum(log_L))
}
