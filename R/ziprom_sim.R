#' Simulate Data from Zero-Inflated Promotion Cure-Rate Models
#'
#' Generates synthetic survival datasets under the zero-adjusted promotion cure-rate model
#' with Gamma frailty and flexible activation schemes (Algorithm 1 of Assuncao et al., 2026).
#'
#' @param n Sample size (integer > 0). Default is 250.
#' @param beta1 Numeric vector of coefficients for zero-inflation linear predictor \eqn{\boldsymbol{\beta}_1}.
#' @param beta2 Numeric vector of coefficients for cure fraction linear predictor \eqn{\boldsymbol{\beta}_2}.
#' @param beta3 Numeric vector of coefficients for shape linear predictor \eqn{\boldsymbol{\beta}_3}.
#' @param beta4 Numeric vector of coefficients for scale linear predictor \eqn{\boldsymbol{\beta}_4}.
#' @param alpha Numeric scalar for Gamma frailty variance (\eqn{\alpha \ge 0}). Default is 0.25.
#' @param activation Character string for activation scheme: \code{"random"} (default),
#'   \code{"minimum"}, or \code{"maximum"}.
#' @param baseline Character string for baseline activation distribution: \code{"weibull"} (default),
#'   \code{"gompertz"}, or \code{"exponential"}.
#' @param censor_rate Nominal right-censoring rate in (0, 1). Default is 0.20.
#' @param seed Optional integer random seed for reproducibility.
#'
#' @return A data frame containing:
#' \describe{
#'   \item{\code{time}}{Observed survival times (\eqn{\min(s_i, w_i)}).}
#'   \item{\code{status}}{Censoring indicator (1 = event, 0 = right-censored).}
#'   \item{\code{x}}{Simulated binary covariate (\eqn{x_i \sim \text{Bernoulli}(0.5)}).}
#' }
#'
#' @references
#' Assuncao, D., Ramos, P. L., Tomazella, V., Leao, J., Osses, L., and Rodrigues, A. (2026).
#' Zero-Inflated Promotion Cure-Rate Models With Frailty and Flexible Activation Schemes.
#' \emph{Biometrical Journal}, 68:e70162. \doi{10.1002/bimj.70162}
#'
#' @export
#' @examples
#' set.seed(123)
#' dat <- sim_ziprom(n = 100, activation = "random")
#' head(dat)
#' table(dat$time == 0) # Instantaneous zero failures
#' table(dat$status)    # Events vs censored
sim_ziprom <- function(n = 250,
                       beta1 = c(-2.0, 0.5),
                       beta2 = c(-1.0, 0.75),
                       beta3 = c(0.5, 0.5),
                       beta4 = c(1.5, 2.0),
                       alpha = 0.25,
                       activation = c("random", "minimum", "maximum"),
                       baseline = c("weibull", "gompertz", "exponential"),
                       censor_rate = 0.20,
                       seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }

  activation <- match.arg(activation)
  baseline <- match.arg(baseline)

  x_cov <- stats::rbinom(n, size = 1, prob = 0.5)
  X1 <- cbind(1, x_cov)
  X2 <- cbind(1, x_cov)
  X3 <- cbind(1, x_cov)
  X4 <- cbind(1, x_cov)

  if (length(beta1) == 1) beta1 <- c(beta1, 0)
  if (length(beta2) == 1) beta2 <- c(beta2, 0)
  if (length(beta3) == 1) beta3 <- c(beta3, 0)
  if (length(beta4) == 1) beta4 <- c(beta4, 0)

  probs <- compute_p0_p1(X1, X2, beta1, beta2)
  p0 <- probs$p0
  p1 <- probs$p1
  theta <- probs$theta

  shape <- exp(as.vector(X3 %*% beta3))
  scale <- exp(as.vector(X4 %*% beta4))

  s_vec <- numeric(n)
  u_vec <- stats::runif(n, 0, 1)

  for (i in 1:n) {
    if (u_vec[i] <= p0[i]) {
      s_vec[i] <- 0
    } else if (u_vec[i] > (1 - p1[i])) {
      s_vec[i] <- Inf
    } else {
      v_i <- (u_vec[i] - p0[i]) / (1 - p0[i] - p1[i])
      s_vec[i] <- quantile_proper(v_i, theta[i], alpha, shape[i], scale[i], activation, baseline)
    }
  }

  pos_finite <- s_vec[is.finite(s_vec) & s_vec > 0]
  tau_censor <- if (length(pos_finite) > 0) stats::quantile(pos_finite, 1 - censor_rate * 0.5) else 10
  c_vec <- stats::runif(n, 0, as.numeric(tau_censor) * 1.5)

  t_vec <- pmin(s_vec, c_vec)
  delta_vec <- ifelse(s_vec <= c_vec, 1, 0)
  delta_vec[s_vec == 0] <- 1

  data.frame(
    time = t_vec,
    status = delta_vec,
    x = x_cov
  )
}

#' Random Sample Generator for Zero-Inflated Promotion Models
#'
#' Alias and wrapper around \code{\link{sim_ziprom}} for generating synthetic survival data.
#'
#' @param n Sample size (integer > 0). Default is 250.
#' @param ... Additional arguments passed to \code{\link{sim_ziprom}}.
#'
#' @return A data frame of simulated observations.
#' @export
#' @examples
#' set.seed(42)
#' d <- rzicr_promotion(n = 50, activation = "minimum")
#' head(d)
rzicr_promotion <- function(n = 250, ...) {
  sim_ziprom(n = n, ...)
}
