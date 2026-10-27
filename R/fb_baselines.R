#' Baseline Hazard and Survival Evaluation
#'
#' Evaluates baseline hazard, cumulative hazard, log hazard, and baseline survival
#' for ten parametric baseline distributions.
#'
#' @param t Non-negative numeric vector of survival times.
#' @param baseline Character specifying baseline distribution: \code{"weibull"},
#'   \code{"exponential"}, \code{"genexponential"}, \code{"exppower"},
#'   \code{"lfr"}, \code{"genloglogistic"}, \code{"genweibull"},
#'   \code{"genrayleigh"}, \code{"weightedexp"}, or \code{"extendedweibull"}.
#' @param params Numeric vector of baseline parameters.
#'
#' @return A list with elements \code{h0} (hazard), \code{H0} (cumulative hazard),
#'   \code{log_h0} (log hazard), and \code{S0} (survival function).
#' @export
#' @keywords internal
eval_baseline <- function(t, baseline = c("weibull", "exponential", "genexponential",
                                          "exppower", "lfr", "genloglogistic",
                                          "genweibull", "genrayleigh", "weightedexp",
                                          "extendedweibull"), params) {
  baseline <- match.arg(baseline)
  t <- pmax(as.numeric(t), 1e-12)
  eps <- 1e-15
  
  if (baseline == "weibull") {
    lambda <- max(params[1], 1e-8)
    rho <- max(params[2], 1e-8)
    H0 <- lambda * (t^rho)
    log_h0 <- log(lambda) + log(rho) + (rho - 1) * log(t)
    h0 <- exp(log_h0)
    S0 <- exp(-H0)
  } else if (baseline == "exponential") {
    lambda <- max(params[1], 1e-8)
    H0 <- lambda * t
    log_h0 <- rep(log(lambda), length(t))
    h0 <- exp(log_h0)
    S0 <- exp(-H0)
  } else if (baseline == "genexponential") {
    lambda <- max(params[1], 1e-8)
    alpha <- max(params[2], 1e-8)
    e_lt <- exp(-lambda * t)
    one_minus_e <- pmax(1 - e_lt, eps)
    F0 <- one_minus_e^alpha
    S0 <- pmax(1 - F0, eps)
    H0 <- -log(S0)
    log_h0 <- log(alpha) + log(lambda) - lambda * t + (alpha - 1) * log(one_minus_e) - log(S0)
    h0 <- exp(log_h0)
  } else if (baseline == "exppower") {
    lambda <- max(params[1], 1e-8)
    alpha <- max(params[2], 1e-8)
    lt_a <- (lambda * t)^alpha
    lt_a_safe <- pmin(lt_a, 700)
    H0 <- exp(lt_a_safe) - 1
    log_h0 <- log(alpha) + alpha * log(lambda) + (alpha - 1) * log(t) + lt_a_safe
    h0 <- exp(log_h0)
    S0 <- exp(-H0)
  } else if (baseline == "lfr") {
    a <- max(params[1], 1e-8)
    b <- max(params[2], 1e-8)
    h0 <- a + b * t
    log_h0 <- log(pmax(h0, eps))
    H0 <- a * t + 0.5 * b * (t^2)
    S0 <- exp(-H0)
  } else if (baseline == "genloglogistic") {
    lambda <- max(params[1], 1e-8)
    alpha <- max(params[2], 1e-8)
    gamma <- max(params[3], 1e-8)
    lt_a <- (lambda * t)^alpha
    one_plus_lta <- 1 + lt_a
    S0 <- pmax(one_plus_lta^(-gamma), eps)
    H0 <- gamma * log(one_plus_lta)
    log_h0 <- log(gamma) + log(alpha) + alpha * log(lambda) + (alpha - 1) * log(t) - log(one_plus_lta)
    h0 <- exp(log_h0)
  } else if (baseline == "genweibull") {
    lambda <- max(params[1], 1e-8)
    alpha <- max(params[2], 1e-8)
    gamma <- max(params[3], 1e-8)
    term <- pmax(1 - lambda * gamma * (t^alpha), eps)
    S0 <- term^(1 / gamma)
    H0 <- -(1 / gamma) * log(term)
    log_h0 <- log(lambda) + log(alpha) + (alpha - 1) * log(t) - log(term)
    h0 <- exp(log_h0)
  } else if (baseline == "genrayleigh") {
    lambda <- max(params[1], 1e-8)
    alpha <- max(params[2], 1e-8)
    lt2 <- (lambda * t)^2
    v <- pmax(1 - exp(-lt2), eps)
    S0 <- pmax(1 - v^alpha, eps)
    H0 <- -log(S0)
    log_h0 <- log(2 * alpha) + 2 * log(lambda) + log(t) - lt2 + (alpha - 1) * log(v) - log(S0)
    h0 <- exp(log_h0)
  } else if (baseline == "weightedexp") {
    lambda <- max(params[1], 1e-8)
    alpha <- max(params[2], 1e-8)
    e1 <- exp(-lambda * t)
    e2 <- exp(-(alpha + 1) * lambda * t)
    S0 <- pmax(((alpha + 1) / alpha) * e1 - (1 / alpha) * e2, eps)
    f0 <- pmax(((alpha + 1) / alpha) * lambda * e1 * (1 - exp(-alpha * lambda * t)), eps)
    H0 <- -log(S0)
    h0 <- f0 / S0
    log_h0 <- log(pmax(h0, eps))
  } else if (baseline == "extendedweibull") {
    lambda <- max(params[1], 1e-8)
    rho <- max(params[2], 1e-8)
    lt_r <- lambda * (t^rho)
    lt_r_safe <- pmin(lt_r, 700)
    H0 <- exp(lt_r_safe) - 1
    log_h0 <- log(lambda) + log(rho) + (rho - 1) * log(t) + lt_r_safe
    h0 <- exp(log_h0)
    S0 <- exp(-H0)
  }
  
  H0 <- pmax(H0, 0)
  S0 <- pmin(pmax(S0, 0), 1)
  
  list(h0 = h0, H0 = H0, log_h0 = log_h0, S0 = S0)
}
