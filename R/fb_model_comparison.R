#' Deviance Information Criterion (DIC)
#'
#' Calculates DIC for a Bayesian frailty model.
#'
#' @param object An object of class \code{"frailty_bayes"}.
#' @return A list with elements \code{DIC}, \code{D_bar}, \code{D_hat}, and \code{p_D}.
#' @export
#' @examples
#' \donttest{
#' data(kidney_frailty)
#' fit <- frailty_bayes(Surv(time1, time2, status1, status2) ~ age,
#'                      data = kidney_frailty, n_iter = 200, burn_in = 50)
#' dic(fit)
#' }
dic <- function(object) {
  if (!inherits(object, "frailty_bayes")) stop("Object must be of class 'frailty_bayes'.")
  list(DIC = object$diagnostics$DIC,
       D_bar = object$diagnostics$D_bar,
       D_hat = object$diagnostics$D_hat,
       p_D = object$diagnostics$p_D)
}

#' Expected Akaike Information Criterion (EAIC)
#' @param object An object of class \code{"frailty_bayes"}.
#' @return Numeric EAIC value.
#' @export
eaic <- function(object) {
  if (!inherits(object, "frailty_bayes")) stop("Object must be of class 'frailty_bayes'.")
  object$diagnostics$EAIC
}

#' Expected Bayesian Information Criterion (EBIC)
#' @param object An object of class \code{"frailty_bayes"}.
#' @return Numeric EBIC value.
#' @export
ebic <- function(object) {
  if (!inherits(object, "frailty_bayes")) stop("Object must be of class 'frailty_bayes'.")
  object$diagnostics$EBIC
}

#' Conditional Predictive Ordinate (CPO) and Log Pseudo Marginal Likelihood (LPML)
#' @param object An object of class \code{"frailty_bayes"}.
#' @return A list containing vector \code{CPO} and scalar \code{LPML}.
#' @export
cpo <- function(object) {
  if (!inherits(object, "frailty_bayes")) stop("Object must be of class 'frailty_bayes'.")
  list(CPO = object$diagnostics$CPO, LPML = object$diagnostics$LPML)
}

#' Log Pseudo Marginal Likelihood (LPML)
#' @param object An object of class \code{"frailty_bayes"}.
#' @return Numeric LPML scalar value.
#' @export
lpml <- function(object) {
  if (!inherits(object, "frailty_bayes")) stop("Object must be of class 'frailty_bayes'.")
  object$diagnostics$LPML
}

#' Bayes Factor Comparison Between Two Models
#'
#' Computes the Bayes factor comparing Model 1 against Model 2 using
#' the ratio of Pseudo Marginal Likelihoods (LPML) or harmonic mean marginal likelihood.
#'
#' @param model1 Fitted object of class \code{"frailty_bayes"}.
#' @param model2 Fitted object of class \code{"frailty_bayes"}.
#' @return A list containing \code{BF12} (Bayes Factor), \code{log_BF12},
#'   and verbal interpretation according to Jeffreys' scale.
#' @export
#' @examples
#' \donttest{
#' data(kidney_frailty)
#' m1 <- frailty_bayes(Surv(time1, time2, status1, status2) ~ age,
#'                     data = kidney_frailty, frailty = "gamma", n_iter = 200, burn_in = 50)
#' m2 <- frailty_bayes(Surv(time1, time2, status1, status2) ~ age,
#'                     data = kidney_frailty, frailty = "invgauss", n_iter = 200, burn_in = 50)
#' bayes_factor(m1, m2)
#' }
bayes_factor <- function(model1, model2) {
  if (!inherits(model1, "frailty_bayes") || !inherits(model2, "frailty_bayes")) {
    stop("Both model1 and model2 must be of class 'frailty_bayes'.")
  }
  
  lpml1 <- model1$diagnostics$LPML
  lpml2 <- model2$diagnostics$LPML
  
  log_bf <- lpml1 - lpml2
  bf12 <- exp(pmin(pmax(log_bf, -700), 700))
  
  interp <- if (log_bf > 2.3) {
    "Decisive evidence in favor of Model 1"
  } else if (log_bf > 1.1) {
    "Strong evidence in favor of Model 1"
  } else if (log_bf > 0.5) {
    "Substantial evidence in favor of Model 1"
  } else if (log_bf > -0.5) {
    "Barely worth mentioning (Inconclusive)"
  } else if (log_bf > -1.1) {
    "Substantial evidence in favor of Model 2"
  } else if (log_bf > -2.3) {
    "Strong evidence in favor of Model 2"
  } else {
    "Decisive evidence in favor of Model 2"
  }
  
  list(BF12 = bf12, log_BF12 = log_bf, interpretation = interp,
       LPML_model1 = lpml1, LPML_model2 = lpml2)
}

#' Model Comparison Table for Frailty Models
#'
#' Creates a unified comparison table ranking multiple fitted frailty models
#' by DIC, EAIC, EBIC, LPML, AIC, and BIC.
#'
#' @param ... One or more fitted objects of class \code{"frailty_bayes"} or \code{"frailty_mle"}.
#' @param names Optional character vector of model names.
#'
#' @return A data frame summarizing model comparison criteria and ranks.
#' @export
#' @examples
#' \donttest{
#' data(kidney_frailty)
#' m1 <- frailty_bayes(Surv(time1, time2, status1, status2) ~ age,
#'                     data = kidney_frailty, frailty = "gamma", n_iter = 200, burn_in = 50)
#' m2 <- frailty_bayes(Surv(time1, time2, status1, status2) ~ age,
#'                     data = kidney_frailty, frailty = "invgauss", n_iter = 200, burn_in = 50)
#' frailty_compare(m1, m2, names = c("Gamma", "InvGauss"))
#' }
frailty_compare <- function(..., names = NULL) {
  models <- list(...)
  n_models <- length(models)
  if (n_models == 0) stop("No models supplied to compare.")
  
  if (is.null(names)) {
    names <- paste0("Model_", 1:n_models)
  }
  
  rows <- vector("list", n_models)
  
  for (i in 1:n_models) {
    m <- models[[i]]
    if (inherits(m, "frailty_bayes")) {
      rows[[i]] <- data.frame(
        Model = names[i],
        Engine = "Bayesian",
        Baseline = m$baseline,
        Frailty = m$frailty,
        Type = m$type,
        Params = m$diagnostics$n_params,
        LogLik = round(m$diagnostics$LogLik, 2),
        DIC = round(m$diagnostics$DIC, 2),
        p_D = round(m$diagnostics$p_D, 2),
        EAIC = round(m$diagnostics$EAIC, 2),
        EBIC = round(m$diagnostics$EBIC, 2),
        LPML = round(m$diagnostics$LPML, 2),
        stringsAsFactors = FALSE
      )
    } else if (inherits(m, "frailty_mle")) {
      rows[[i]] <- data.frame(
        Model = names[i],
        Engine = "MLE",
        Baseline = m$baseline,
        Frailty = m$frailty,
        Type = m$type,
        Params = m$n_params,
        LogLik = round(m$logLik, 2),
        DIC = NA,
        p_D = NA,
        EAIC = round(m$AIC, 2),
        EBIC = round(m$BIC, 2),
        LPML = NA,
        stringsAsFactors = FALSE
      )
    }
  }
  
  res <- do.call(rbind, rows)
  if (all(!is.na(res$DIC))) {
    res$Rank_DIC <- rank(res$DIC)
  }
  res
}
