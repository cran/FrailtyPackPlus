#' Baseline and Frailty Distributions for Additive Shared Frailty Models
#'
#' Evaluates survival functions, hazard rates, cumulative hazard functions,
#' frailty probability densities, Laplace transforms, and derivatives for
#' Generalized Weibull (GW), Generalized Log-Logistic Type-II (GLL2),
#' Generalized Lindley (GL), Weighted Lindley (WL), Gamma, and Inverse Gaussian.
#'
#' @param t Numeric vector of positive survival times.
#' @param delta Positive scale parameter (> 0).
#' @param zeta Positive shape/exponent parameter (> 0).
#' @param xi Positive power parameter (> 0).
#' @param w,v Numeric vector of positive frailty values (> 0).
#' @param s Numeric vector of non-negative Laplace transform arguments (>= 0).
#' @param theta Positive frailty parameter (> 0).
#' @param mu Positive frailty parameter with 0 < mu < 1 + theta.
#' @param s_param Positive shape parameter for Weighted Lindley frailty (> 0).
#'
#' @return A numeric vector of calculated distribution values or Laplace transform evaluations.
#'
#' @name distributions
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
#' Mudholkar, G. S., Srivastava, D. K., and Freimer, M. (1995). The exponentiated
#' Weibull family: A reanalysis of the bus-motor-failure data.
#' \emph{Technometrics}, 37(4), 436-445. \doi{10.1080/00401706.1995.10484376}
#'
#' @examples
#' # Generalized Weibull baseline survival and hazard
#' s_gw(c(10, 20, 30), delta = 0.05, zeta = 2.0, xi = 0.7)
#' h_gw(c(10, 20, 30), delta = 0.05, zeta = 2.0, xi = 0.7)
#'
#' # Generalized Lindley Laplace transform and derivatives
#' laplace_genlindley(c(0.1, 0.5), theta = 2.5, mu = 1.5)
#' laplace_deriv1_genlindley(c(0.1, 0.5), theta = 2.5, mu = 1.5)
#'
#' # Weighted Lindley Laplace transform
#' laplace_wlindley(c(0.1, 0.5), s_param = 1.5)
NULL

#' @rdname distributions
#' @export
s_gw <- function(t, delta, zeta, xi) {
  t <- pmax(t, 0)
  u <- -expm1(-delta * (t^xi))
  u <- pmin(pmax(u, 0), 1 - 1e-15)
  1 - (u^zeta)
}

#' @rdname distributions
#' @export
cumh_gw <- function(t, delta, zeta, xi) {
  surv <- s_gw(t, delta, zeta, xi)
  surv <- pmax(surv, 1e-15)
  -log(surv)
}

#' @rdname distributions
#' @export
h_gw <- function(t, delta, zeta, xi) {
  t <- pmax(t, 1e-10)
  t_xi <- t^xi
  exp_term <- exp(-delta * t_xi)
  u <- -expm1(-delta * t_xi)
  u <- pmin(pmax(u, 1e-15), 1 - 1e-15)
  surv <- pmax(1 - (u^zeta), 1e-15)
  num <- zeta * delta * xi * (t^(xi - 1)) * exp_term * (u^(zeta - 1))
  num / surv
}

#' @rdname distributions
#' @export
s_gll2 <- function(t, delta, zeta, xi) {
  t <- pmax(t, 0)
  (1 + delta * (t^xi))^(-zeta)
}

#' @rdname distributions
#' @export
cumh_gll2 <- function(t, delta, zeta, xi) {
  t <- pmax(t, 0)
  zeta * log1p(delta * (t^xi))
}

#' @rdname distributions
#' @export
h_gll2 <- function(t, delta, zeta, xi) {
  t <- pmax(t, 1e-10)
  (zeta * delta * xi * (t^(xi - 1))) / (1 + delta * (t^xi))
}

#' @rdname distributions
#' @export
d_genlindley <- function(w, theta, mu) {
  if (theta <= 0 || mu <= 0 || mu >= 1 + theta) return(rep(0, length(w)))
  eta <- theta * (1 + theta - mu)
  w_pos <- pmax(w, 1e-15)
  term1 <- (theta^(mu + 1) * (w_pos^(mu - 1))) / gamma(mu)
  term2 <- (theta^eta * (w_pos^(eta - 1))) / gamma(eta)
  (1 / (1 + theta)) * (term1 + term2) * exp(-theta * w_pos)
}

#' @rdname distributions
#' @export
laplace_genlindley <- function(s, theta, mu) {
  s <- pmax(s, 0)
  eta <- theta * (1 + theta - mu)
  denom <- theta + s
  c1 <- (theta^(mu + 1)) / (1 + theta)
  c2 <- (theta^eta) / (1 + theta)
  c1 * (denom^(-mu)) + c2 * (denom^(-eta))
}

#' @rdname distributions
#' @export
laplace_deriv1_genlindley <- function(s, theta, mu) {
  s <- pmax(s, 0)
  eta <- theta * (1 + theta - mu)
  denom <- theta + s
  c1 <- (theta^(mu + 1)) / (1 + theta)
  c2 <- (theta^eta) / (1 + theta)
  -c1 * mu * (denom^(-(mu + 1))) - c2 * eta * (denom^(-(eta + 1)))
}

#' @rdname distributions
#' @export
laplace_deriv2_genlindley <- function(s, theta, mu) {
  s <- pmax(s, 0)
  eta <- theta * (1 + theta - mu)
  denom <- theta + s
  c1 <- (theta^(mu + 1)) / (1 + theta)
  c2 <- (theta^eta) / (1 + theta)
  c1 * mu * (mu + 1) * (denom^(-(mu + 2))) + c2 * eta * (eta + 1) * (denom^(-(eta + 2)))
}

#' @rdname distributions
#' @export
var_genlindley <- function(theta, mu) {
  if (theta <= 0 || mu <= 0 || mu >= 1 + theta) return(NA_real_)
  num <- theta^4 - (theta^3) * mu + 3 * (theta^2) * (1 + theta) - 4 * (theta^2) * mu + 3 * theta * mu * (mu - 1) + mu^2
  den <- theta * ((1 + theta)^2)
  num / den
}

#' @rdname distributions
#' @export
d_wlindley <- function(v, s_param) {
  if (s_param <= 0) return(rep(0, length(v)))
  c_val <- sqrt(s_param * (s_param + 1))
  v_pos <- pmax(v, 1e-15)
  coef <- (c_val^(s_param + 1)) / ((c_val + s_param) * gamma(s_param))
  coef * (v_pos^(s_param - 1)) * (1 + v_pos) * exp(-c_val * v_pos)
}

#' @rdname distributions
#' @export
laplace_wlindley <- function(s, s_param) {
  s <- pmax(s, 0)
  c_val <- sqrt(s_param * (s_param + 1))
  c0 <- (c_val^(s_param + 1)) / (c_val + s_param)
  denom <- c_val + s
  c0 * (denom^(-s_param) + s_param * (denom^(-(s_param + 1))))
}

#' @rdname distributions
#' @export
laplace_deriv1_wlindley <- function(s, s_param) {
  s <- pmax(s, 0)
  c_val <- sqrt(s_param * (s_param + 1))
  c0 <- (c_val^(s_param + 1)) / (c_val + s_param)
  denom <- c_val + s
  -c0 * (s_param * (denom^(-(s_param + 1))) + s_param * (s_param + 1) * (denom^(-(s_param + 2))))
}

#' @rdname distributions
#' @export
laplace_deriv2_wlindley <- function(s, s_param) {
  s <- pmax(s, 0)
  c_val <- sqrt(s_param * (s_param + 1))
  c0 <- (c_val^(s_param + 1)) / (c_val + s_param)
  denom <- c_val + s
  c0 * (s_param * (s_param + 1) * (denom^(-(s_param + 2))) + s_param * (s_param + 1) * (s_param + 2) * (denom^(-(s_param + 3))))
}

#' @rdname distributions
#' @export
var_wlindley <- function(s_param) {
  if (s_param <= 0) return(NA_real_)
  2 / (s_param + sqrt(s_param * (s_param + 1)))
}

#' @rdname distributions
#' @export
laplace_gamma <- function(s, theta) {
  s <- pmax(s, 0)
  (1 + theta * s)^(-1 / theta)
}

#' @rdname distributions
#' @export
laplace_deriv1_gamma <- function(s, theta) {
  s <- pmax(s, 0)
  -(1 + theta * s)^(-1 / theta - 1)
}

#' @rdname distributions
#' @export
laplace_deriv2_gamma <- function(s, theta) {
  s <- pmax(s, 0)
  (1 + theta) * ((1 + theta * s)^(-1 / theta - 2))
}

#' @rdname distributions
#' @export
var_gamma <- function(theta) {
  if (theta <= 0) return(NA_real_)
  theta
}

#' @rdname distributions
#' @export
laplace_invgaussian <- function(s, theta) {
  s <- pmax(s, 0)
  exp((1 / theta) * (1 - sqrt(1 + 2 * theta * s)))
}

#' @rdname distributions
#' @export
laplace_deriv1_invgaussian <- function(s, theta) {
  s <- pmax(s, 0)
  sqrt_term <- sqrt(1 + 2 * theta * s)
  - (1 / sqrt_term) * exp((1 / theta) * (1 - sqrt_term))
}

#' @rdname distributions
#' @export
laplace_deriv2_invgaussian <- function(s, theta) {
  s <- pmax(s, 0)
  sqrt_term <- sqrt(1 + 2 * theta * s)
  lt <- exp((1 / theta) * (1 - sqrt_term))
  ( (1 + theta * sqrt_term) / (sqrt_term^3) ) * lt
}

#' @rdname distributions
#' @export
var_invgaussian <- function(theta) {
  if (theta <= 0) return(NA_real_)
  theta
}
