# Model Comparison and Bayes Factor Functions

#' Compare Competing Shared RHR Frailty Models
#'
#' Tabulates and compares multiple fitted Bayesian or MLE shared frailty models
#' using information criteria (AIC, BIC, AICc, HQIC, DIC, LPML), Kendall's tau, and frailty variance.
#'
#' @details
#' Compares goodness-of-fit and penalized likelihood criteria across competing models.
#' Models with lower AIC, BIC, AICc, HQIC, and DIC, and higher LPML (closer to zero) represent superior fits.
#'
#' @param ... Fitted objects of class \code{"rhr_bayes_fit"} or \code{"rhr_mle_fit"}.
#' @param models_list Optional list of fitted model objects.
#'
#' @return A data frame with one row per model containing:
#' \describe{
#'   \item{\code{Model}}{Name or label of the model.}
#'   \item{\code{Frailty}}{Frailty distribution.}
#'   \item{\code{Baseline}}{Baseline distribution.}
#'   \item{\code{Method}}{Estimation methodology (\code{"Bayes"} or \code{"MLE"}).}
#'   \item{\code{LogLik}}{Log-likelihood at parameter point estimates.}
#'   \item{\code{AIC}}{Akaike Information Criterion.}
#'   \item{\code{BIC}}{Bayesian Information Criterion.}
#'   \item{\code{AICc}}{Corrected Akaike Information Criterion.}
#'   \item{\code{HQIC}}{Hannan-Quinn Information Criterion.}
#'   \item{\code{DIC}}{Deviance Information Criterion (for Bayesian fits).}
#'   \item{\code{LPML}}{Log Pseudo Marginal Likelihood (for Bayesian fits).}
#'   \item{\code{Kendall_Tau}}{Estimated Kendall's tau rank correlation.}
#'   \item{\code{Frailty_Var}}{Estimated theoretical variance of the frailty variable.}
#' }
#' @export
#'
#' @examples
#' data(austrian_twin)
#' fit1 <- fit_rhr_frailty_bayes(cbind(t1, t2, d1, d2) ~ age, data = austrian_twin[1:30, ],
#'                               frailty = "WL", baseline = "GE", n_iter = 40, burn_in = 10)
#' fit2 <- fit_rhr_frailty_bayes(cbind(t1, t2, d1, d2) ~ age, data = austrian_twin[1:30, ],
#'                               frailty = "WL", baseline = "EG", n_iter = 40, burn_in = 10)
#' compare_rhr_models(WL_GE = fit1, WL_EG = fit2)
compare_rhr_models <- function(..., models_list = NULL) {
  dots <- list(...)
  if (!is.null(models_list)) dots <- c(dots, models_list)
  if (length(dots) == 0) stop("No models provided for comparison.")
  
  mod_names <- names(dots)
  if (is.null(mod_names) || any(mod_names == "")) {
    mod_names <- paste0("Model_", seq_along(dots))
  }
  
  rows <- vector("list", length(dots))
  
  for (i in seq_along(dots)) {
    m <- dots[[i]]
    is_bayes <- inherits(m, "rhr_bayes_fit")
    is_mle <- inherits(m, "rhr_mle_fit")
    
    if (!is_bayes && !is_mle) {
      warning("Model ", mod_names[i], " is neither rhr_bayes_fit nor rhr_mle_fit. Skipped.")
      next
    }
    
    rows[[i]] <- data.frame(
      Model = mod_names[i],
      Frailty = m$model_info$frailty,
      Baseline = m$model_info$baseline,
      Method = if (is_bayes) "Bayes" else "MLE",
      LogLik = if (is_bayes) round(-0.5 * m$criteria$dev_at_mean, 2) else round(m$loglik, 2),
      AIC = round(m$criteria$AIC, 2),
      BIC = round(m$criteria$BIC, 2),
      AICc = round(m$criteria$AICc, 2),
      HQIC = round(m$criteria$HQIC, 2),
      DIC = if (is_bayes) round(m$criteria$DIC, 2) else NA_real_,
      LPML = if (is_bayes) round(m$criteria$LPML, 2) else NA_real_,
      Kendall_Tau = round(m$kendall_tau, 4),
      Frailty_Var = round(m$frailty_variance, 4),
      stringsAsFactors = FALSE
    )
  }
  
  res_df <- do.call(rbind, rows)
  rownames(res_df) <- NULL
  return(res_df)
}

#' Bayes Factor Calculation for Model Comparison
#'
#' Computes the Bayes factor comparing two Bayesian shared RHR frailty models
#' based on log pseudo-marginal likelihoods (LPML) or marginal deviances:
#' \deqn{BF_{12} = \exp(LPML_1 - LPML_2)}
#'
#' @details
#' Evaluates the strength of evidence between two Bayesian models using the scale of Kass and Raftery (1995)
#' based on \eqn{2 \ln(BF_{12})}.
#'
#' @param model1 Fitted object of class \code{"rhr_bayes_fit"}.
#' @param model2 Fitted object of class \code{"rhr_bayes_fit"}.
#'
#' @return A list containing the following components:
#' \describe{
#'   \item{\code{log_BF12}}{Natural logarithm of the Bayes factor \eqn{\ln(BF_{12})}.}
#'   \item{\code{BF12}}{Bayes factor value \eqn{BF_{12}}.}
#'   \item{\code{two_log_BF12}}{Twice the natural logarithm of the Bayes factor \eqn{2 \ln(BF_{12})}.}
#'   \item{\code{LPML_Model1}}{Log pseudo-marginal likelihood of Model 1.}
#'   \item{\code{LPML_Model2}}{Log pseudo-marginal likelihood of Model 2.}
#'   \item{\code{Evidence}}{Verbal interpretation of evidence according to the Kass and Raftery (1995) scale.}
#' }
#' @export
#'
#' @examples
#' data(austrian_twin)
#' fit1 <- fit_rhr_frailty_bayes(cbind(t1, t2, d1, d2) ~ age, data = austrian_twin[1:30, ],
#'                               frailty = "WL", baseline = "GE", n_iter = 40, burn_in = 10)
#' fit2 <- fit_rhr_frailty_bayes(cbind(t1, t2, d1, d2) ~ age, data = austrian_twin[1:30, ],
#'                               frailty = "WL", baseline = "EG", n_iter = 40, burn_in = 10)
#' bayes_factor_rhr(fit1, fit2)
bayes_factor_rhr <- function(model1, model2) {
  if (!inherits(model1, "rhr_bayes_fit") || !inherits(model2, "rhr_bayes_fit")) {
    stop("Both models must be of class 'rhr_bayes_fit'.")
  }
  
  lpml1 <- model1$criteria$LPML
  lpml2 <- model2$criteria$LPML
  
  log_bf <- lpml1 - lpml2
  bf12 <- exp(pmin(500, pmax(-500, log_bf)))
  
  two_log_bf <- 2 * log_bf
  evidence <- if (two_log_bf > 10) {
    "Very strong evidence for Model 1"
  } else if (two_log_bf > 6) {
    "Strong evidence for Model 1"
  } else if (two_log_bf > 2) {
    "Positive evidence for Model 1"
  } else if (two_log_bf > 0) {
    "Weak / barely noteworthy evidence for Model 1"
  } else if (two_log_bf > -2) {
    "Weak / barely noteworthy evidence for Model 2"
  } else if (two_log_bf > -6) {
    "Positive evidence for Model 2"
  } else if (two_log_bf > -10) {
    "Strong evidence for Model 2"
  } else {
    "Very strong evidence for Model 2"
  }
  
  list(
    log_BF12 = log_bf,
    BF12 = bf12,
    two_log_BF12 = two_log_bf,
    LPML_Model1 = lpml1,
    LPML_Model2 = lpml2,
    Evidence = evidence
  )
}
