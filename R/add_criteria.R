#' Model Selection Criteria and Bayes Factor Comparisons
#'
#' Computes AIC, BIC, Corrected AIC (AICc), Hannan-Quinn Criterion (HQIC),
#' Deviance Information Criterion (DIC), and Bayes Factors for additive
#' shared frailty models.
#'
#' @param loglik Scalar maximized or posterior mean log-likelihood value.
#' @param k Number of estimated model parameters.
#' @param n Sample size (number of subjects/pairs).
#' @param deviance_samples Optional numeric vector of posterior deviances (-2 * logLik) for DIC calculation.
#' @param deviance_at_mean Optional scalar deviance evaluated at the posterior mean of parameters.
#' @param fit1,fit2 Fitted model objects of class \code{"add_frailty"}.
#'
#' @return For \code{calc_model_criteria}, a named list of information criteria.
#'   For \code{calc_bayes_factor}, a scalar Bayes factor comparing \code{fit1} against \code{fit2}.
#'
#' @name criteria
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
#' # Calculate information criteria
#' calc_model_criteria(loglik = -330.5, k = 13, n = 38)
NULL

#' @rdname criteria
#' @export
calc_model_criteria <- function(loglik, k, n, deviance_samples = NULL, deviance_at_mean = NULL) {
  aic <- 2 * k - 2 * loglik
  bic <- k * log(n) - 2 * loglik
  aicc <- if (n > k + 1) aic + (2 * k * (k + 1)) / (n - k - 1) else aic
  hqic <- 2 * k * log(log(pmax(n, 3))) - 2 * loglik

  dic <- NA_real_
  if (!is.null(deviance_samples) && !is.null(deviance_at_mean)) {
    valid_dev <- deviance_samples[is.finite(deviance_samples) & !is.nan(deviance_samples)]
    if (length(valid_dev) > 0 && is.finite(deviance_at_mean)) {
      d_bar <- mean(valid_dev)
      dic <- 2 * d_bar - deviance_at_mean
    }
  }

  list(
    LogLik = loglik,
    k = k,
    n = n,
    AIC = aic,
    BIC = bic,
    AICc = aicc,
    HQIC = hqic,
    DIC = dic
  )
}

#' @rdname criteria
#' @export
calc_bayes_factor <- function(fit1, fit2) {
  get_marg_loglik <- function(fit) {
    if (!is.null(fit$marginal_loglik)) {
      return(fit$marginal_loglik)
    }
    if (!is.null(fit$deviance_samples)) {
      # Newton-Raftery harmonic mean estimator
      ll_samples <- -0.5 * fit$deviance_samples
      max_ll <- max(ll_samples)
      - (log(mean(exp(- (ll_samples - max_ll)))) - max_ll)
    } else if (!is.null(fit$loglik) && !is.null(fit$k) && !is.null(fit$n)) {
      # BIC approximation to marginal log-likelihood
      fit$loglik - 0.5 * fit$k * log(fit$n)
    } else {
      stop("Unable to extract or compute marginal likelihood from fitted object.")
    }
  }

  ml1 <- get_marg_loglik(fit1)
  ml2 <- get_marg_loglik(fit2)

  exp(pmin(pmax(ml1 - ml2, -50), 50))
}
