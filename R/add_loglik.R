#' Bivariate Survival Log-Likelihood for Additive Shared Frailty Models
#'
#' Evaluates the exact log-likelihood function under random right censoring
#' for additive shared frailty regression models with Generalized Lindley,
#' Weighted Lindley, Gamma, or Inverse Gaussian frailty, and Generalized Weibull
#' or Generalized Log-Logistic Type-II baseline distributions.
#'
#' @param params Numeric vector of all model parameters:
#'   \code{c(zeta1, delta1, xi1, zeta2, delta2, xi2, frailty_params, beta)}.
#' @param t1 Numeric vector of first survival times (positive).
#' @param delta1 Numeric vector of first censoring indicators (1 = event, 0 = censored).
#' @param t2 Numeric vector of second survival times (positive).
#' @param delta2 Numeric vector of second censoring indicators (1 = event, 0 = censored).
#' @param X Covariate design matrix with \code{nrow(X) == length(t1)}.
#' @param frailty_type Character string specifying the frailty distribution:
#'   \code{"genlindley"} (Generalized Lindley), \code{"wlindley"} (Weighted Lindley),
#'   \code{"gamma"} (Gamma), or \code{"invgaussian"} (Inverse Gaussian).
#' @param baseline_dist Character string specifying baseline distribution:
#'   \code{"gw"} (Generalized Weibull) or \code{"gll2"} (Generalized Log-Logistic Type-II).
#'
#' @return A scalar numeric value representing the total log-likelihood.
#'
#' @export
#' @references
#' Pandey, A., Hanagal, D. D., and Tyagi, S. (2022). Generalised Lindley shared
#' additive frailty regression model for bivariate survival data.
#' \emph{Statistics in Transition new series}, 23(4), 161-176.
#' \doi{10.2478/stattrans-2022-0048}
#'
#' Tyagi, S., Pandey, A., Hanagal, D. D., and Chesneau, C. (2024). A study on
#' comparisons of additive regression frailty models to counter heterogeneity:
#' Bayesian strategies and case study. \emph{Communications in Statistics -
#' Simulation and Computation}. \doi{10.1080/03610918.2024.2391874}
#'
#' @examples
#' data(kidney_bivariate)
#' X <- as.matrix(kidney_bivariate[, c("age", "sex", "gn", "an", "pkd")])
#' pars <- c(zeta1 = 2.0, delta1 = 0.05, xi1 = 0.7,
#'           zeta2 = 2.7, delta2 = 0.06, xi2 = 0.7,
#'           theta = 2.3, mu = 1.4,
#'           beta1 = -0.1, beta2 = -9.0, beta3 = 2.4, beta4 = 1.6, beta5 = -50.0)
#' loglik_add_frailty(pars, kidney_bivariate$t1, kidney_bivariate$delta1,
#'                    kidney_bivariate$t2, kidney_bivariate$delta2,
#'                    X, frailty_type = "genlindley", baseline_dist = "gw")
loglik_add_frailty <- function(params, t1, delta1, t2, delta2, X,
                               frailty_type = c("genlindley", "wlindley", "gamma", "invgaussian"),
                               baseline_dist = c("gw", "gll2")) {
  frailty_type <- match.arg(frailty_type)
  baseline_dist <- match.arg(baseline_dist)

  zeta1 <- params[1]
  delta1_par <- params[2]
  xi1 <- params[3]
  zeta2 <- params[4]
  delta2_par <- params[5]
  xi2 <- params[6]

  if (any(c(zeta1, delta1_par, xi1, zeta2, delta2_par, xi2) <= 0)) {
    return(-1e12)
  }

  if (frailty_type == "genlindley") {
    theta <- params[7]
    mu <- params[8]
    if (theta <= 0 || mu <= 0 || mu >= 1 + theta) return(-1e12)
    beta <- params[9:length(params)]
  } else if (frailty_type == "wlindley") {
    s_param <- params[7]
    if (s_param <= 0) return(-1e12)
    beta <- params[8:length(params)]
  } else if (frailty_type %in% c("gamma", "invgaussian")) {
    theta <- params[7]
    if (theta <= 0) return(-1e12)
    beta <- params[8:length(params)]
  }

  if (is.null(X) || ncol(X) == 0) {
    rho <- rep(1, length(t1))
  } else {
    lin_pred <- as.vector(as.matrix(X) %*% beta)
    lin_pred <- pmin(pmax(lin_pred, -30), 30)
    rho <- exp(lin_pred)
  }

  s_arg <- (t1 + t2) * rho

  if (baseline_dist == "gw") {
    Phi01 <- cumh_gw(t1, delta1_par, zeta1, xi1)
    Phi02 <- cumh_gw(t2, delta2_par, zeta2, xi2)
    phi01 <- h_gw(t1, delta1_par, zeta1, xi1)
    phi02 <- h_gw(t2, delta2_par, zeta2, xi2)
  } else {
    Phi01 <- cumh_gll2(t1, delta1_par, zeta1, xi1)
    Phi02 <- cumh_gll2(t2, delta2_par, zeta2, xi2)
    phi01 <- h_gll2(t1, delta1_par, zeta1, xi1)
    phi02 <- h_gll2(t2, delta2_par, zeta2, xi2)
  }

  if (frailty_type == "genlindley") {
    L0 <- laplace_genlindley(s_arg, theta, mu)
    L1 <- laplace_deriv1_genlindley(s_arg, theta, mu)
    L2 <- laplace_deriv2_genlindley(s_arg, theta, mu)
  } else if (frailty_type == "wlindley") {
    L0 <- laplace_wlindley(s_arg, s_param)
    L1 <- laplace_deriv1_wlindley(s_arg, s_param)
    L2 <- laplace_deriv2_wlindley(s_arg, s_param)
  } else if (frailty_type == "gamma") {
    L0 <- laplace_gamma(s_arg, theta)
    L1 <- laplace_deriv1_gamma(s_arg, theta)
    L2 <- laplace_deriv2_gamma(s_arg, theta)
  } else {
    L0 <- laplace_invgaussian(s_arg, theta)
    L1 <- laplace_deriv1_invgaussian(s_arg, theta)
    L2 <- laplace_deriv2_invgaussian(s_arg, theta)
  }

  term <- numeric(length(t1))

  c11 <- (delta1 == 1 & delta2 == 1)
  c10 <- (delta1 == 1 & delta2 == 0)
  c01 <- (delta1 == 0 & delta2 == 1)
  c00 <- (delta1 == 0 & delta2 == 0)

  if (any(c11)) {
    term[c11] <- phi01[c11] * phi02[c11] * L0[c11] -
                 rho[c11] * (phi01[c11] + phi02[c11]) * L1[c11] +
                 (rho[c11]^2) * L2[c11]
  }
  if (any(c10)) {
    term[c10] <- phi01[c10] * L0[c10] - rho[c10] * L1[c10]
  }
  if (any(c01)) {
    term[c01] <- phi02[c01] * L0[c01] - rho[c01] * L1[c01]
  }
  if (any(c00)) {
    term[c00] <- L0[c00]
  }

  term <- pmax(term, 1e-300)
  ind_loglik <- -(Phi01 + Phi02) + log(term)

  if (any(is.nan(ind_loglik) | is.na(ind_loglik))) {
    return(-1e12)
  }

  sum(ind_loglik)
}
