#' Simulate Frailty Survival Data or One-Shot Device Lifetimes
#'
#' @description
#' Simulates synthetic survival datasets under shared, additive, or multiplicative
#' frailty models, or simulates correlated component lifetimes for one-shot devices
#' under shared gamma frailty.
#'
#' @param n An integer specifying sample size (number of subjects, clusters, or devices). Default is 100.
#' @param baseline Baseline hazard distribution name: \code{"weibull"} (default),
#'   \code{"exponential"}, \code{"exppower"}, etc.
#' @param frailty Frailty distribution name: \code{"gamma"} (default), etc.
#' @param type Model type: \code{"multiplicative_bivariate"} (default),
#'   \code{"additive_bivariate"}, \code{"multiplicative_univariate"}, or \code{"additive_univariate"}.
#' @param beta Regression coefficients vector or frailty variance parameter (for one-shot devices).
#' @param theta Frailty parameter(s).
#' @param baseline_params Baseline hazard parameter vector.
#' @param censor_rate Proportion of right-censored observations (default 0.2).
#' @param lambda Optional positive numeric vector of component failure rates \eqn{\lambda_m}
#'   (for one-shot device testing simulation). When supplied, one-shot device lifetimes are simulated.
#' @param alpha Optional positive numeric vector of Weibull shape parameters (for one-shot devices).
#' @param seed Optional integer for reproducibility.
#'
#' @return If \code{lambda} is supplied, an \eqn{n \times M} numeric matrix of simulated
#'   one-shot component lifetimes. Otherwise, a data frame containing simulated survival times,
#'   status indicators, covariates, and true frailty values. Missing values are not generated.
#' @export
#' @examples
#' # Simulate survival data
#' sim_dat <- sim_frailty(n = 50, baseline = "weibull", frailty = "gamma")
#' head(sim_dat)
#'
#' # Simulate one-shot device lifetimes
#' lam <- c(0.0086, 0.0067, 0.0052, 0.0024)
#' sim_times <- sim_frailty(n = 50, lambda = lam, beta = 0.3, seed = 123)
#' head(sim_times)
sim_frailty <- function(n = 100, baseline = "weibull", frailty = "gamma",
                        type = "multiplicative_bivariate", beta = c(0.5, -0.5),
                        theta = 1.0, baseline_params = c(0.01, 1.2),
                        censor_rate = 0.2, lambda = NULL, alpha = NULL, seed = NULL) {
  if (!is.null(lambda)) {
    return(sim_oneshot_frailty_lifetimes(n = n, lambda = lambda, beta = beta, alpha = alpha, seed = seed))
  }
  if (!is.null(seed)) set.seed(seed)
  fb_sim_frailty(n = n, baseline = baseline, frailty = frailty, type = type,
                 beta = beta, theta = theta, baseline_params = baseline_params,
                 censor_rate = censor_rate)
}

#' Simulate One-Shot Device Lifetimes under Shared Gamma Frailty
#'
#' @description
#' Convenience alias to generate correlated component lifetimes for one-shot devices under gamma frailty.
#'
#' @param n An integer specifying the number of devices to generate.
#' @param lambda A positive numeric vector of failure rates \eqn{\lambda_m} for each component.
#' @param beta A positive numeric scalar representing the frailty variance (\eqn{\beta > 0}).
#' @param alpha An optional positive numeric vector of Weibull shape parameters.
#' @param seed An optional integer for reproducibility.
#'
#' @return An \eqn{n \times M} numeric matrix of simulated component failure times. Missing values are not permitted.
#' @export
sim_oneshot_frailty <- function(n, lambda, beta, alpha = NULL, seed = NULL) {
  sim_oneshot_frailty_lifetimes(n = n, lambda = lambda, beta = beta, alpha = alpha, seed = seed)
}
