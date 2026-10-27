#' Cumulative Baseline Hazard Function
#'
#' Computes the cumulative baseline hazard \eqn{H_0(t)} for Weibull, Gompertz,
#' or Exponential baseline distributions.
#'
#' @param t Numeric vector of positive survival times.
#' @param shape Numeric vector or scalar representing baseline shape parameter(s) (> 0).
#' @param scale Numeric vector or scalar representing baseline scale parameter(s) (> 0).
#' @param baseline Character string indicating the baseline hazard distribution:
#'   \code{"weibull"} (default), \code{"gompertz"}, or \code{"exponential"}.
#'
#' @return Numeric vector of cumulative baseline hazard values.
#' @export
base_cumhaz <- function(t, shape, scale, baseline = c("weibull", "gompertz", "exponential")) {
  baseline <- match.arg(baseline)
  t <- pmax(t, 1e-12)
  shape <- pmax(shape, 1e-8)
  scale <- pmax(scale, 1e-8)

  if (baseline == "weibull") {
    (t / scale)^shape
  } else if (baseline == "gompertz") {
    (shape / scale) * (exp(pmin(scale * t, 500)) - 1)
  } else if (baseline == "exponential") {
    t / scale
  }
}

#' Baseline Hazard Function
#'
#' Computes the baseline hazard rate \eqn{h_0(t)} for Weibull, Gompertz,
#' or Exponential baseline distributions.
#'
#' @param t Numeric vector of positive survival times.
#' @param shape Numeric vector or scalar representing baseline shape parameter(s) (> 0).
#' @param scale Numeric vector or scalar representing baseline scale parameter(s) (> 0).
#' @param baseline Character string indicating the baseline hazard distribution:
#'   \code{"weibull"} (default), \code{"gompertz"}, or \code{"exponential"}.
#'
#' @return Numeric vector of baseline hazard rate values.
#' @export
base_haz <- function(t, shape, scale, baseline = c("weibull", "gompertz", "exponential")) {
  baseline <- match.arg(baseline)
  t <- pmax(t, 1e-12)
  shape <- pmax(shape, 1e-8)
  scale <- pmax(scale, 1e-8)

  if (baseline == "weibull") {
    (shape / scale) * (t / scale)^(shape - 1)
  } else if (baseline == "gompertz") {
    shape * exp(pmin(scale * t, 500))
  } else if (baseline == "exponential") {
    rep(1 / scale, length(t))
  }
}

#' Compute Multinomial Logit Probabilities for Zero-Inflation and Cure Fraction
#'
#' Evaluates the probabilities of zero-failure \eqn{p_0}, cure fraction \eqn{p_1},
#' susceptible fraction \eqn{1 - p_0 - p_1}, and the latent Poisson parameter \eqn{\theta = -\log(p_1)}.
#'
#' @param X1 Design matrix for zero-inflation linear predictor.
#' @param X2 Design matrix for cure fraction linear predictor.
#' @param beta1 Numeric vector of coefficients for zero-inflation.
#' @param beta2 Numeric vector of coefficients for cure fraction.
#'
#' @return A list containing \code{p0}, \code{p1}, \code{theta}, and \code{susceptible}.
#' @export
compute_p0_p1 <- function(X1, X2, beta1, beta2) {
  eta1 <- as.vector(X1 %*% beta1)
  eta2 <- as.vector(X2 %*% beta2)

  max_eta <- pmax(0, eta1, eta2)
  exp0 <- exp(-max_eta)
  exp1 <- exp(eta1 - max_eta)
  exp2 <- exp(eta2 - max_eta)
  denom <- exp0 + exp1 + exp2

  p0 <- exp1 / denom
  p1 <- exp2 / denom

  p0 <- pmin(pmax(p0, 1e-12), 1 - 1e-6)
  p1 <- pmin(pmax(p1, 1e-12), 1 - p0 - 1e-6)
  theta <- -log(pmax(p1, 1e-12))
  susceptible <- pmax(1 - p0 - p1, 1e-12)

  list(p0 = p0, p1 = p1, theta = theta, susceptible = susceptible)
}

#' Proper Conditional Survival Function for Susceptible Group
#'
#' Computes \eqn{S_p^*(t) = P(T > t \mid M \ge 1)} under Minimum, Maximum,
#' and Random activation schemes with Gamma frailty.
#'
#' @param t Numeric vector of positive survival times.
#' @param theta Numeric vector or scalar representing the latent Poisson parameter (> 0).
#' @param alpha Numeric scalar representing Gamma frailty variance (>= 0).
#' @param shape Numeric vector or scalar representing baseline shape parameter(s) (> 0).
#' @param scale Numeric vector or scalar representing baseline scale parameter(s) (> 0).
#' @param activation Character string indicating the activation scheme:
#'   \code{"random"} (default), \code{"minimum"}, or \code{"maximum"}.
#' @param baseline Character string indicating the baseline hazard distribution:
#'   \code{"weibull"} (default), \code{"gompertz"}, or \code{"exponential"}.
#'
#' @return Numeric vector of proper conditional survival probabilities.
#' @export
surv_proper <- function(t, theta, alpha, shape, scale,
                        activation = c("random", "minimum", "maximum"),
                        baseline = c("weibull", "gompertz", "exponential")) {
  activation <- match.arg(activation)
  baseline <- match.arg(baseline)

  H0 <- base_cumhaz(t, shape, scale, baseline)
  theta <- pmax(theta, 1e-8)

  if (alpha > 1e-8) {
    u <- (1 + alpha * H0)^(-1 / alpha)
  } else {
    u <- exp(-pmin(H0, 500))
  }
  u <- pmin(pmax(u, 1e-300), 1)

  if (activation == "random") {
    res <- u
  } else if (activation == "minimum") {
    num <- exp(pmin(-theta * (1 - u), 0)) - exp(-theta)
    denom <- pmax(1 - exp(-theta), 1e-12)
    res <- num / denom
  } else if (activation == "maximum") {
    num <- 1 - exp(-pmin(theta * u, 500))
    denom <- pmax(1 - exp(-theta), 1e-12)
    res <- num / denom
  }

  pmin(pmax(res, 0), 1)
}

#' Proper Conditional Density Function for Susceptible Group
#'
#' Computes \eqn{f_p^*(t)} for \eqn{t > 0} under Minimum, Maximum,
#' and Random activation schemes with Gamma frailty.
#'
#' @param t Numeric vector of positive survival times.
#' @param theta Numeric vector or scalar representing the latent Poisson parameter (> 0).
#' @param alpha Numeric scalar representing Gamma frailty variance (>= 0).
#' @param shape Numeric vector or scalar representing baseline shape parameter(s) (> 0).
#' @param scale Numeric vector or scalar representing baseline scale parameter(s) (> 0).
#' @param activation Character string indicating the activation scheme:
#'   \code{"random"} (default), \code{"minimum"}, or \code{"maximum"}.
#' @param baseline Character string indicating the baseline hazard distribution:
#'   \code{"weibull"} (default), \code{"gompertz"}, or \code{"exponential"}.
#'
#' @return Numeric vector of proper conditional density values.
#' @export
dens_proper <- function(t, theta, alpha, shape, scale,
                        activation = c("random", "minimum", "maximum"),
                        baseline = c("weibull", "gompertz", "exponential")) {
  activation <- match.arg(activation)
  baseline <- match.arg(baseline)

  H0 <- base_cumhaz(t, shape, scale, baseline)
  h0 <- base_haz(t, shape, scale, baseline)
  theta <- pmax(theta, 1e-8)

  if (alpha > 1e-8) {
    u <- (1 + alpha * H0)^(-1 / alpha)
    du <- (1 + alpha * H0)^(-(alpha + 1) / alpha) * h0
  } else {
    u <- exp(-pmin(H0, 500))
    du <- h0 * u
  }
  u <- pmin(pmax(u, 1e-300), 1)
  du <- pmax(du, 1e-300)

  if (activation == "random") {
    res <- du
  } else if (activation == "minimum") {
    denom <- pmax(1 - exp(-theta), 1e-12)
    res <- (theta * du * exp(pmin(-theta * (1 - u), 0))) / denom
  } else if (activation == "maximum") {
    denom <- pmax(1 - exp(-theta), 1e-12)
    res <- (theta * du * exp(-pmin(theta * u, 500))) / denom
  }

  pmax(res, 1e-300)
}

#' Population Survival Function
#'
#' Computes the unconditional population survival function \eqn{S_p(t)}.
#'
#' @param t Numeric vector of non-negative times.
#' @param p0 Numeric vector or scalar of zero-failure probabilities.
#' @param p1 Numeric vector or scalar of cure fractions.
#' @param alpha Numeric scalar representing Gamma frailty variance (>= 0).
#' @param shape Numeric vector or scalar for shape parameter.
#' @param scale Numeric vector or scalar for scale parameter.
#' @param activation Character string for activation scheme (\code{"random"}, \code{"minimum"}, or \code{"maximum"}).
#' @param baseline Character string for baseline hazard distribution (\code{"weibull"}, \code{"gompertz"}, or \code{"exponential"}).
#'
#' @return Numeric vector of population survival probabilities.
#' @export
#' @examples
#' t_seq <- seq(0, 10, by = 1)
#' s_vals <- ziprom_core_surv(t_seq, p0 = 0.1, p1 = 0.4, alpha = 0.25,
#'                           shape = 1.2, scale = 2.5, activation = "random")
#' print(s_vals)
ziprom_core_surv <- function(t, p0, p1, alpha = 0, shape = 1, scale = 1,
                             activation = c("random", "minimum", "maximum"),
                             baseline = c("weibull", "gompertz", "exponential")) {
  activation <- match.arg(activation)
  baseline <- match.arg(baseline)

  theta <- -log(pmax(p1, 1e-12))
  susceptible <- pmax(1 - p0 - p1, 0)

  s_star <- ifelse(t == 0, 1, surv_proper(t, theta, alpha, shape, scale, activation, baseline))
  res <- ifelse(t == 0, 1 - p0, p1 + susceptible * s_star)
  pmin(pmax(res, 0), 1)
}

#' Population Density Function
#'
#' Computes the population density function \eqn{f_p(t)} for positive survival times \eqn{t > 0}.
#'
#' @param t Numeric vector of positive times (\eqn{t > 0}).
#' @param p0 Numeric vector or scalar of zero-failure probabilities.
#' @param p1 Numeric vector or scalar of cure fractions.
#' @param alpha Numeric scalar representing Gamma frailty variance (>= 0).
#' @param shape Numeric vector or scalar for shape parameter.
#' @param scale Numeric vector or scalar for scale parameter.
#' @param activation Character string for activation scheme (\code{"random"}, \code{"minimum"}, or \code{"maximum"}).
#' @param baseline Character string for baseline hazard distribution (\code{"weibull"}, \code{"gompertz"}, or \code{"exponential"}).
#'
#' @return Numeric vector of population density values for \eqn{t > 0}.
#' @export
#' @examples
#' t_seq <- seq(0.1, 10, by = 1)
#' d_vals <- ziprom_core_dens(t_seq, p0 = 0.1, p1 = 0.4, alpha = 0.25,
#'                           shape = 1.2, scale = 2.5, activation = "random")
#' print(d_vals)
ziprom_core_dens <- function(t, p0, p1, alpha = 0, shape = 1, scale = 1,
                             activation = c("random", "minimum", "maximum"),
                             baseline = c("weibull", "gompertz", "exponential")) {
  activation <- match.arg(activation)
  baseline <- match.arg(baseline)

  theta <- -log(pmax(p1, 1e-12))
  susceptible <- pmax(1 - p0 - p1, 0)

  f_star <- dens_proper(t, theta, alpha, shape, scale, activation, baseline)
  susceptible * f_star
}

#' Quantile Function for Susceptible Group
#'
#' Computes the quantile survival time \eqn{t} such that \eqn{S_p^*(t) = 1 - v}
#' for uniform probabilities \eqn{v \in (0, 1)}.
#'
#' @param v Numeric vector of probabilities in (0, 1).
#' @param theta Numeric vector or scalar representing the latent Poisson parameter (> 0).
#' @param alpha Numeric scalar representing Gamma frailty variance (>= 0).
#' @param shape Numeric vector or scalar representing baseline shape parameter(s) (> 0).
#' @param scale Numeric vector or scalar representing baseline scale parameter(s) (> 0).
#' @param activation Character string indicating the activation scheme:
#'   \code{"random"} (default), \code{"minimum"}, or \code{"maximum"}.
#' @param baseline Character string indicating the baseline hazard distribution:
#'   \code{"weibull"} (default), \code{"gompertz"}, or \code{"exponential"}.
#'
#' @return Numeric vector of quantile survival times for susceptible individuals.
#' @export
quantile_proper <- function(v, theta, alpha, shape, scale,
                            activation = c("random", "minimum", "maximum"),
                            baseline = c("weibull", "gompertz", "exponential")) {
  activation <- match.arg(activation)
  baseline <- match.arg(baseline)

  v <- pmin(pmax(v, 1e-10), 1 - 1e-10)
  theta <- pmax(theta, 1e-8)

  if (activation == "random") {
    w <- 1 - v
  } else if (activation == "minimum") {
    inner <- exp(-theta) + (1 - v) * (1 - exp(-theta))
    w <- 1 + (1 / theta) * log(pmax(inner, 1e-300))
  } else if (activation == "maximum") {
    inner <- 1 - (1 - v) * (1 - exp(-theta))
    w <- -(1 / theta) * log(pmax(inner, 1e-300))
  }
  w <- pmin(pmax(w, 1e-12), 1)

  if (alpha > 1e-8) {
    H0 <- (w^(-alpha) - 1) / alpha
  } else {
    H0 <- -log(w)
  }
  H0 <- pmax(H0, 0)

  if (baseline == "weibull") {
    scale * (H0)^(1 / shape)
  } else if (baseline == "gompertz") {
    (1 / scale) * log(pmax(1 + (scale / shape) * H0, 1e-12))
  } else if (baseline == "exponential") {
    scale * H0
  }
}
