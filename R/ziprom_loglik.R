#' Log-Likelihood Function for Zero-Inflated Promotion Cure-Rate Models
#'
#' Computes the exact log-likelihood for the zero-adjusted promotion cure-rate model
#' with Gamma frailty and activation schemes under random right censoring.
#'
#' @param par Numeric vector of parameters: \code{c(beta1, beta2, beta3, beta4, log_alpha)}
#'   where \code{beta1} is for zero inflation, \code{beta2} is for cure fraction,
#'   \code{beta3} is for baseline shape, \code{beta4} is for baseline scale,
#'   and \code{log_alpha} is log(frailty variance) if \code{frailty = TRUE}.
#' @param time Numeric vector of non-negative observed survival times.
#' @param status Numeric vector of event indicators (1 = event, 0 = censored).
#' @param X1 Design matrix for zero-inflation probability \eqn{p_0}.
#' @param X2 Design matrix for cure probability \eqn{p_1}.
#' @param X3 Design matrix for baseline shape parameter \eqn{\alpha_1}.
#' @param X4 Design matrix for baseline scale parameter \eqn{\alpha_2}.
#' @param activation Character string: \code{"random"}, \code{"minimum"}, or \code{"maximum"}.
#' @param baseline Character string: \code{"weibull"}, \code{"gompertz"}, or \code{"exponential"}.
#' @param frailty Logical indicating whether Gamma frailty is included (\code{TRUE} or \code{FALSE}).
#'
#' @return Scalar value of the negative log-likelihood (for minimization in \code{optim}),
#'   or the positive log-likelihood if called with appropriate sign.
#'
#' @references
#' Assuncao, D., Ramos, P. L., Tomazella, V., Leao, J., Osses, L., and Rodrigues, A. (2026).
#' Zero-Inflated Promotion Cure-Rate Models With Frailty and Flexible Activation Schemes.
#' \emph{Biometrical Journal}, 68:e70162. \doi{10.1002/bimj.70162}
#'
#' @export
#' @examples
#' set.seed(123)
#' dat <- sim_ziprom(n = 50, activation = "random")
#' X <- matrix(1, nrow = nrow(dat), ncol = 1)
#' par_init <- c(-0.5, 1.0, 0.0, 1.0, log(0.25))
#' nll <- ziprom_loglik(par_init, dat$time, dat$status, X, X, X, X,
#'                      activation = "random", baseline = "weibull", frailty = TRUE)
#' print(nll)
ziprom_loglik <- function(par, time, status, X1, X2, X3, X4,
                          activation = c("random", "minimum", "maximum"),
                          baseline = c("weibull", "gompertz", "exponential"),
                          frailty = TRUE) {
  activation <- match.arg(activation)
  baseline <- match.arg(baseline)

  p1_len <- ncol(X1)
  p2_len <- ncol(X2)
  p3_len <- ncol(X3)
  p4_len <- ncol(X4)

  idx1 <- 1:p1_len
  idx2 <- (p1_len + 1):(p1_len + p2_len)
  idx3 <- (p1_len + p2_len + 1):(p1_len + p2_len + p3_len)
  idx4 <- (p1_len + p2_len + p3_len + 1):(p1_len + p2_len + p3_len + p4_len)

  beta1 <- par[idx1]
  beta2 <- par[idx2]
  beta3 <- par[idx3]
  beta4 <- par[idx4]

  if (isTRUE(frailty)) {
    log_alpha <- par[length(par)]
    alpha <- exp(pmin(pmax(log_alpha, -20), 20))
  } else {
    alpha <- 0
  }

  probs <- compute_p0_p1(X1, X2, beta1, beta2)
  p0 <- probs$p0
  p1 <- probs$p1
  theta <- probs$theta
  susceptible <- probs$susceptible

  shape <- exp(pmin(pmax(as.vector(X3 %*% beta3), -10), 10))
  scale <- exp(pmin(pmax(as.vector(X4 %*% beta4), -10), 10))

  # Subgroups
  is_zero <- (time == 0)
  is_event <- (!is_zero & status == 1)
  is_censored <- (!is_zero & status == 0)

  ll_total <- 0

  # 1. Zeros: log(p0)
  if (any(is_zero)) {
    ll_total <- ll_total + sum(log(pmax(p0[is_zero], 1e-300)))
  }

  # 2. Events (t > 0, delta = 1): log(1 - p0 - p1) + log(f_p^*(t))
  if (any(is_event)) {
    t_ev <- time[is_event]
    th_ev <- theta[is_event]
    sh_ev <- shape[is_event]
    sc_ev <- scale[is_event]
    susc_ev <- susceptible[is_event]

    f_star <- dens_proper(t_ev, th_ev, alpha, sh_ev, sc_ev, activation, baseline)
    ll_total <- ll_total + sum(log(pmax(susc_ev, 1e-300)) + log(pmax(f_star, 1e-300)))
  }

  # 3. Censored (t > 0, delta = 0): log(p1 + (1 - p0 - p1) * S_p^*(t))
  if (any(is_censored)) {
    t_ce <- time[is_censored]
    th_ce <- theta[is_censored]
    sh_ce <- shape[is_censored]
    sc_ce <- scale[is_censored]
    p1_ce <- p1[is_censored]
    susc_ce <- susceptible[is_censored]

    s_star <- surv_proper(t_ce, th_ce, alpha, sh_ce, sc_ce, activation, baseline)
    surv_pop <- p1_ce + susc_ce * s_star
    ll_total <- ll_total + sum(log(pmax(surv_pop, 1e-300)))
  }

  if (is.na(ll_total) || is.nan(ll_total) || is.infinite(ll_total)) {
    return(1e10)
  }

  -ll_total
}
