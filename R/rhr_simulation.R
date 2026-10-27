# Simulation of Bivariate Left-Censored Data under Shared RHR Frailty Models

#' Simulate Bivariate Left-Censored Survival Data under Shared RHR Frailty
#'
#' Generates synthetic bivariate left-censored survival data with covariates and shared frailty
#' based on reversed hazard rates.
#'
#' @details
#' Simulates the latent frailty variable \eqn{W_j} from the specified frailty distribution,
#' then generates conditionally independent lifetimes \eqn{T_{1j}, T_{2j}} via probability integral transforms on conditional baseline reversed hazard distributions:
#' \deqn{F(t \mid W_j, X_j) = [F_0(t)]^{W_j \exp(X_j^T \\beta)}}
#' Left-censoring is applied such that if \eqn{T_{ij} \le C_{ij}}, the recorded time is \eqn{C_{ij}} and \eqn{\delta_{ij} = 0}; otherwise \eqn{\delta_{ij} = 1}.
#' Missing values (\code{NA} or \code{NaN}) are not generated.
#'
#' @param n Sample size (number of pairs). Default is \code{100}.
#' @param frailty Character string specifying frailty distribution. Default is \code{"WL"}.
#' @param baseline Character string specifying baseline distribution. Default is \code{"GE"}.
#' @param theta Numeric parameter for frailty distribution. Default is \code{3.4}.
#' @param mu Optional secondary parameter for \code{"GL"}.
#' @param theta2 Optional second parameter for \code{"GIG"}.
#' @param params1 Numeric named vector of baseline parameters for first component. Default is \code{c(zeta = 12.4, delta = 0.11)}.
#' @param params2 Numeric named vector of baseline parameters for second component. Default is \code{c(zeta = 11.8, delta = 0.11)}.
#' @param beta Numeric vector of regression coefficients. Default is \code{c(0.0075, -0.11, 0.13)}.
#' @param covars Optional numeric matrix of covariates of dimension \code{n x p}. If \code{NULL}, covariates are generated automatically.
#' @param censor_time1 Numeric scalar or vector of left-censoring thresholds for T1. Default is \code{11}.
#' @param censor_time2 Numeric scalar or vector of left-censoring thresholds for T2. Default is \code{11}.
#' @param seed Optional integer seed for reproducibility. Default is \code{NULL}.
#'
#' @return A data frame containing simulated columns:
#' \describe{
#'   \item{\code{t1}}{Observed or left-censored lifetime for component 1.}
#'   \item{\code{t2}}{Observed or left-censored lifetime for component 2.}
#'   \item{\code{d1}}{Censoring indicator for component 1 (\code{1} = observed, \code{0} = left-censored).}
#'   \item{\code{d2}}{Censoring indicator for component 2 (\code{1} = observed, \code{0} = left-censored).}
#'   \item{covariates}{Columns corresponding to the generated or supplied covariates.}
#'   \item{\code{frailty_w}}{The unobserved latent frailty value \eqn{W_j}.}
#' }
#' @export
#'
#' @examples
#' sim_data <- sim_biv_rhr_frailty(n = 50, frailty = "WL", baseline = "GE", seed = 123)
#' head(sim_data)
sim_biv_rhr_frailty <- function(n = 100, frailty = c("WL", "GL", "GIG", "Gamma", "IG", "Lindley", "None"),
                                baseline = c("GE", "EG", "GIR", "MIW", "GLL-I", "GLL-II", "GIE", "Weibull", "Exponential", "Log-Logistic"),
                                theta = 3.4, mu = NULL, theta2 = NULL,
                                params1 = c(zeta = 12.4, delta = 0.11),
                                params2 = c(zeta = 11.8, delta = 0.11),
                                beta = c(0.0075, -0.11, 0.13),
                                covars = NULL, censor_time1 = 11, censor_time2 = 11, seed = NULL) {
  frailty <- match.arg(frailty)
  baseline <- match.arg(baseline)
  
  if (!is.null(seed)) set.seed(seed)
  
  if (is.null(covars)) {
    age <- stats::runif(n, 18, 55)
    k2 <- stats::rbinom(n, 1, 0.45)
    k3 <- stats::rbinom(n, 1, 0.42)
    covars <- cbind(age = age, k2 = k2, k3 = k3)
  }
  
  p_cov <- ncol(covars)
  if (length(beta) != p_cov) {
    beta <- rep(0.1, p_cov)
  }
  
  eta <- exp(as.numeric(covars %*% beta))
  
  if (frailty == "WL") {
    C <- sqrt(theta * (theta + 1))
    p_mix <- C / (C + theta)
    comp <- stats::rbinom(n, 1, 1 - p_mix)
    w <- stats::rgamma(n, shape = theta + comp, rate = C)
  } else if (frailty == "GL") {
    if (is.null(mu)) mu <- 0.5 * (1 + theta)
    eta_gl <- theta * (1 + theta - mu)
    comp <- stats::rbinom(n, 1, 1 / (1 + theta))
    w <- numeric(n)
    w[comp == 1] <- stats::rgamma(sum(comp == 1), shape = mu, rate = theta)
    w[comp == 0] <- stats::rgamma(sum(comp == 0), shape = eta_gl, rate = theta)
  } else if (frailty == "GIG") {
    th2 <- if (is.null(theta2)) theta else theta2
    w1 <- theta / (theta + th2)
    comp <- stats::rbinom(n, 1, 1 - w1)
    rig <- function(nn, mu_val, lambda_val) {
      y <- stats::rnorm(nn)^2
      x <- mu_val + (mu_val^2 * y) / (2 * lambda_val) - (mu_val / (2 * lambda_val)) * sqrt(4 * mu_val * lambda_val * y + mu_val^2 * y^2)
      z <- stats::runif(nn)
      ifelse(z <= mu_val / (mu_val + x), x, mu_val^2 / x)
    }
    w <- numeric(n)
    w[comp == 0] <- rig(sum(comp == 0), 1.0, theta)
    w[comp == 1] <- rig(sum(comp == 1), 1.0, th2)
  } else if (frailty == "Gamma") {
    w <- stats::rgamma(n, shape = theta, rate = theta)
  } else if (frailty == "IG") {
    rig <- function(nn, mu_val, lambda_val) {
      y <- stats::rnorm(nn)^2
      x <- mu_val + (mu_val^2 * y) / (2 * lambda_val) - (mu_val / (2 * lambda_val)) * sqrt(4 * mu_val * lambda_val * y + mu_val^2 * y^2)
      z <- stats::runif(nn)
      ifelse(z <= mu_val / (mu_val + x), x, mu_val^2 / x)
    }
    w <- rig(n, 1.0, theta)
  } else if (frailty == "Lindley") {
    comp <- stats::rbinom(n, 1, 1 / (theta + 1))
    w <- stats::rgamma(n, shape = 1 + comp, rate = theta)
  } else {
    w <- rep(1.0, n)
  }
  
  u1 <- stats::runif(n)
  u2 <- stats::runif(n)
  
  p1_adj <- u1^(1 / (w * eta))
  p2_adj <- u2^(1 / (w * eta))
  
  t1_raw <- q_rhr_baseline(p1_adj, baseline, params1)
  t2_raw <- q_rhr_baseline(p2_adj, baseline, params2)
  
  c1_vec <- if (length(censor_time1) == 1) rep(censor_time1, n) else censor_time1
  c2_vec <- if (length(censor_time2) == 1) rep(censor_time2, n) else censor_time2
  
  d1 <- ifelse(t1_raw >= c1_vec, 1L, 0L)
  t1 <- ifelse(d1 == 1L, t1_raw, c1_vec)
  
  d2 <- ifelse(t2_raw >= c2_vec, 1L, 0L)
  t2 <- ifelse(d2 == 1L, t2_raw, c2_vec)
  
  sim_df <- data.frame(
    t1 = round(t1, 2),
    t2 = round(t2, 2),
    d1 = d1,
    d2 = d2,
    covars,
    frailty_w = round(w, 4)
  )
  
  return(sim_df)
}
