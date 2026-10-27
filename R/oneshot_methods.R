`%||%` <- function(a, b) if (!is.null(a)) a else b

#' Print Method for One-Shot Frailty Model Fits
#'
#' @description
#' Prints a concise summary of the fitted \code{"oneshot_kofm"} model object.
#'
#' @param x A fitted model object of class \code{"oneshot_kofm"}.
#' @param digits An integer specifying the number of significant digits (default 4).
#' @param ... Additional arguments.
#'
#' @return No return value, called for printing side effects.
#'
#' @export
print.oneshot_kofm <- function(x, digits = 4, ...) {
  cat("\n===================================================================\n")
  cat(" Inference for k-Out-of-M One-Shot Devices with Gamma Frailty\n")
  cat("===================================================================\n\n")
  cat("Call:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  
  cat("Model Specifications:\n")
  cat("  Components (M):         ", x$M, " (", paste(x$components, collapse = ", "), ")\n", sep = "")
  cat("  Baseline Distribution:  ", toupper(substr(x$baseline, 1, 1)), substr(x$baseline, 2, nchar(x$baseline)), "\n", sep = "")
  cat("  Frailty Formulation:    ", if (x$frailty) "Shared Gamma Frailty" else "Independence Model", "\n", sep = "")
  cat("  Estimation Algorithm:   ", toupper(x$method), "\n", sep = "")
  cat("  Total Sample Size:      ", x$sample_size, " devices across ", nrow(x$fitted_probs), " test groups\n\n", sep = "")
  
  cat("Coefficients:\n")
  print(round(x$coefficients, digits))
  cat("\nLog-Likelihood: ", format(round(x$loglik, digits), nsmall = digits),
      "  (df = ", x$n_params, ")\n", sep = "")
  cat("AIC: ", format(round(x$aic, digits), nsmall = digits),
      "   BIC: ", format(round(x$bic, digits), nsmall = digits), "\n", sep = "")
  
  if (x$frailty && !is.null(x$kendall_tau)) {
    cat("Kendall's tau association: ", format(round(x$kendall_tau$estimate, digits), nsmall = digits),
        " (SE = ", format(round(x$kendall_tau$se, digits), nsmall = digits), ")\n", sep = "")
  }
  cat("\n")
  invisible(x)
}

#' Summary Method for One-Shot Frailty Model Fits
#'
#' @description
#' Computes and formats comprehensive estimation summaries, standard errors,
#' \eqn{z}-values, \eqn{p}-values, and confidence intervals for \code{"oneshot_kofm"} objects.
#'
#' @param object A fitted model object of class \code{"oneshot_kofm"}.
#' @param ci_level A numeric scalar between 0 and 1 specifying the confidence level (default from object).
#' @param ... Additional arguments.
#'
#' @return An S3 object of class \code{"summary.oneshot_kofm"} containing formatted tables
#'   and model evaluation statistics.
#'
#' @export
summary.oneshot_kofm <- function(object, ci_level = object$ci_level, ...) {
  cf <- object$coefficients
  se <- object$se
  z_val <- cf / ifelse(se > 0, se, NA)
  p_val <- 2 * (1 - stats::pnorm(abs(z_val)))
  
  # Wald test for beta is one-sided
  if (object$frailty && "beta" %in% names(cf)) {
    p_val["beta"] <- stats::pnorm(z_val["beta"], lower.tail = FALSE)
  }
  
  z_crit <- stats::qnorm(1 - (1 - ci_level) / 2)
  lci <- cf - z_crit * se
  uci <- cf + z_crit * se
  if (object$frailty && "beta" %in% names(cf)) {
    lci["beta"] <- max(0, lci["beta"])
    uci["beta"] <- min(0.5, max(0, uci["beta"]))
  }
  
  coef_table <- cbind(
    Estimate = cf,
    `Std. Error` = se,
    `z value` = z_val,
    `Pr(>|z|)` = p_val,
    `LCI` = lci,
    `UCI` = uci
  )
  
  res <- list(
    call = object$call,
    coefficients = coef_table,
    loglik = object$loglik,
    n_params = object$n_params,
    aic = object$aic,
    bic = object$bic,
    aicc = object$aicc,
    hqic = object$hqic,
    sample_size = object$sample_size,
    M = object$M,
    components = object$components,
    baseline = object$baseline,
    frailty = object$frailty,
    method = object$method,
    kendall_tau = object$kendall_tau,
    ci_level = ci_level
  )
  class(res) <- "summary.oneshot_kofm"
  res
}

#' Print Summary for One-Shot Frailty Model Fits
#'
#' @description
#' Prints the summary table for class \code{"summary.oneshot_kofm"}.
#'
#' @param x An object of class \code{"summary.oneshot_kofm"}.
#' @param digits An integer specifying the number of significant digits (default 4).
#' @param ... Additional arguments.
#'
#' @return No return value, called for printing side effects.
#'
#' @export
print.summary.oneshot_kofm <- function(x, digits = 4, ...) {
  cat("\n===================================================================\n")
  cat(" Inference for k-Out-of-M One-Shot Devices with Gamma Frailty\n")
  cat("===================================================================\n\n")
  cat("Call:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  
  cat("Model Specifications:\n")
  cat("  Components (M):         ", x$M, " (", paste(x$components, collapse = ", "), ")\n", sep = "")
  cat("  Baseline Distribution:  ", toupper(substr(x$baseline, 1, 1)), substr(x$baseline, 2, nchar(x$baseline)), "\n", sep = "")
  cat("  Frailty Formulation:    ", if (x$frailty) "Shared Gamma Frailty" else "Independence Model", "\n", sep = "")
  cat("  Estimation Algorithm:   ", toupper(x$method), "\n", sep = "")
  cat("  Sample Size:            ", x$sample_size, " devices\n\n", sep = "")
  
  cat("Coefficients & Asymptotic ", x$ci_level * 100, "% Confidence Intervals:\n", sep = "")
  stats::printCoefmat(x$coefficients, digits = digits, cs.ind = 1:2, tst.ind = 3, p.column = 4, has.Pvalue = TRUE, signif.stars = FALSE)
  
  cat("\nModel Selection Statistics:\n")
  cat("  Log-Likelihood: ", format(round(x$loglik, digits), nsmall = digits), " (df = ", x$n_params, ")\n", sep = "")
  cat("  AIC:            ", format(round(x$aic, digits), nsmall = digits), "\n", sep = "")
  cat("  AICc:           ", format(round(x$aicc, digits), nsmall = digits), "\n", sep = "")
  cat("  BIC:            ", format(round(x$bic, digits), nsmall = digits), "\n", sep = "")
  cat("  HQIC:           ", format(round(x$hqic, digits), nsmall = digits), "\n", sep = "")
  
  if (x$frailty && !is.null(x$kendall_tau)) {
    cat("\nDependence Association:\n")
    cat("  Kendall's tau:  ", format(round(x$kendall_tau$estimate, digits), nsmall = digits),
        " (SE = ", format(round(x$kendall_tau$se, digits), nsmall = digits), ", 95% CI: [",
        format(round(x$kendall_tau$confint["LCI"], digits), nsmall = digits), ", ",
        format(round(x$kendall_tau$confint["UCI"], digits), nsmall = digits), "])\n", sep = "")
  }
  cat("\n")
  invisible(x)
}

#' Extract Model Coefficients
#'
#' @description
#' Extracts estimated parameters from a fitted \code{"oneshot_kofm"} object.
#'
#' @param object A fitted model object of class \code{"oneshot_kofm"}.
#' @param ... Additional arguments.
#'
#' @return A named numeric vector of parameter estimates.
#'
#' @export
coef.oneshot_kofm <- function(object, ...) {
  object$coefficients
}

#' Extract Variance-Covariance Matrix
#'
#' @description
#' Extracts the asymptotic variance-covariance matrix of parameter estimators
#' from a fitted \code{"oneshot_kofm"} object.
#'
#' @param object A fitted model object of class \code{"oneshot_kofm"}.
#' @param ... Additional arguments.
#'
#' @return A symmetric numeric variance-covariance matrix.
#'
#' @export
vcov.oneshot_kofm <- function(object, ...) {
  object$vcov
}

#' Confidence Intervals for Model Parameters
#'
#' @description
#' Computes asymptotic confidence intervals for parameters of a \code{"oneshot_kofm"} model.
#'
#' @param object A fitted model object of class \code{"oneshot_kofm"}.
#' @param parm A specification of which parameters are to be given confidence intervals.
#' @param level The confidence level required (default 0.95).
#' @param ... Additional arguments.
#'
#' @return A numeric matrix with columns \code{LCI} and \code{UCI}.
#'
#' @export
confint.oneshot_kofm <- function(object, parm, level = 0.95, ...) {
  cf <- object$coefficients
  se <- object$se
  z_crit <- stats::qnorm(1 - (1 - level) / 2)
  
  lci <- cf - z_crit * se
  uci <- cf + z_crit * se
  if (object$frailty && "beta" %in% names(cf)) {
    lci["beta"] <- max(0, lci["beta"])
    uci["beta"] <- min(0.5, max(0, uci["beta"]))
  }
  
  ci_mat <- cbind(LCI = lci, UCI = uci)
  rownames(ci_mat) <- names(cf)
  
  if (!missing(parm)) {
    ci_mat <- ci_mat[parm, , drop = FALSE]
  }
  ci_mat
}

#' Extract Model Log-Likelihood
#'
#' @description
#' Extracts the maximized log-likelihood value from a fitted \code{"oneshot_kofm"} object.
#'
#' @param object A fitted model object of class \code{"oneshot_kofm"}.
#' @param ... Additional arguments.
#'
#' @return An object of class \code{"logLik"} with degrees of freedom attribute.
#'
#' @export
logLik.oneshot_kofm <- function(object, ...) {
  structure(
    object$loglik,
    df = object$n_params,
    nobs = object$sample_size,
    class = "logLik"
  )
}

#' Extract Fitted Values
#'
#' @description
#' Extracts fitted cell probabilities or expected counts from a \code{"oneshot_kofm"} object.
#'
#' @param object A fitted model object of class \code{"oneshot_kofm"}.
#' @param type A character string: \code{"probs"} (default) or \code{"counts"}.
#' @param ... Additional arguments.
#'
#' @return A numeric matrix of fitted probabilities or expected counts.
#'
#' @export
fitted.oneshot_kofm <- function(object, type = c("probs", "counts"), ...) {
  type <- match.arg(type)
  if (type == "probs") object$fitted_probs else object$fitted_counts
}

#' Extract Pearson Residuals
#'
#' @description
#' Extracts standardized Pearson residuals from a fitted \code{"oneshot_kofm"} object.
#'
#' @param object A fitted model object of class \code{"oneshot_kofm"}.
#' @param type A character string: \code{"pearson"} (default).
#' @param ... Additional arguments.
#'
#' @return A numeric matrix of standardized Pearson residuals.
#'
#' @export
residuals.oneshot_kofm <- function(object, type = c("pearson"), ...) {
  type <- match.arg(type)
  object$residuals
}

#' Predict System Reliability or Hazard Rates
#'
#' @description
#' Generates predictions of \eqn{k}-out-of-\eqn{M} system reliability or component failure rates
#' for new stress levels and time intervals from a fitted \code{"oneshot_kofm"} object.
#'
#' @param object A fitted model object of class \code{"oneshot_kofm"}.
#' @param newdata An optional data frame with new \code{stress} or \code{time} values.
#' @param type A character string: \code{"reliability"} (default), \code{"mttf"}, or \code{"moments"}.
#' @param k An optional integer vector specifying system structures.
#' @param ... Additional arguments passed to underlying prediction functions.
#'
#' @return A data frame of predictions.
#'
#' @export
predict.oneshot_kofm <- function(object, newdata = NULL,
                                 type = c("reliability", "mttf", "moments"),
                                 k = NULL, ...) {
  type <- match.arg(type)
  
  if (is.null(newdata)) {
    newdata <- object$data
  }
  
  if (type == "reliability") {
    stress_val <- if ("stress" %in% colnames(newdata)) newdata$stress else newdata[, 1]
    time_val <- if ("time" %in% colnames(newdata)) newdata$time else 1000
    reliability_k_out_of_m(object, t = time_val, stress = stress_val, k = k)
  } else if (type == "mttf") {
    stress_val <- if ("stress" %in% colnames(newdata)) newdata$stress else newdata[, 1]
    mttf_k_out_of_m(object, stress = stress_val, k = k, ...)
  } else {
    stress_val <- if ("stress" %in% colnames(newdata)) newdata$stress else newdata[, 1]
    component_moments(object, stress = stress_val)
  }
}


# S3 Methods for oneshot_frailty

`%||%` <- function(a, b) if (!is.null(a)) a else b

#' @export
print.oneshot_frailty <- function(x, digits = 4, ...) {
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  cat(sprintf("One-Shot Device Frailty Model (%s Baseline)\n",
              if (x$baseline == "weibull") "Weibull" else "Exponential"))
  cat(sprintf("Dependence Structure: %s\n",
              if (x$frailty) sprintf("Shared Gamma Frailty (beta = %.*f, Kendall's tau = %.*f)",
                                     digits, x$coefficients["beta"], digits, x$kendall_tau$tau)
              else "Independent Components (beta = 0)"))
  cat(sprintf("Components (%d): %s\n", x$M, paste(x$components, collapse = ", ")))
  cat(sprintf("Log-Likelihood: %.*f (df = %d)\n", digits, x$loglik, x$n_params))
  cat(sprintf("AIC: %.*f, BIC: %.*f\n\n", digits, x$aic, digits, x$bic))
  
  cat("Coefficients:\n")
  print(round(x$coefficients, digits))
  invisible(x)
}

#' @export
summary.oneshot_frailty <- function(object, ...) {
  coefs <- object$coefficients
  se <- object$se
  z_vals <- coefs / se
  p_vals <- 2 * stats::pnorm(abs(z_vals), lower.tail = FALSE)
  p_vals <- pmin(pmax(as.numeric(p_vals), 0), 1)
  p_vals[is.na(p_vals) | is.nan(p_vals)] <- 1
  
  ci_low <- object$confint[, "LCI"]
  ci_high <- object$confint[, "UCI"]
  
  coef_table <- cbind(
    Estimate = coefs,
    `Std. Error` = se,
    `z value` = z_vals,
    `Pr(>|z|)` = p_vals
  )
  
  ci_table <- cbind(
    `95% LCI` = ci_low,
    `95% UCI` = ci_high
  )
  
  res <- list(
    call = object$call,
    baseline = object$baseline,
    frailty = object$frailty,
    M = object$M,
    components = object$components,
    coefficients = coef_table,
    confint = ci_table,
    loglik = object$loglik,
    n_params = object$n_params,
    aic = object$aic,
    bic = object$bic,
    aicc = object$aicc,
    hqic = object$hqic,
    sample_size = object$sample_size,
    kendall_tau = object$kendall_tau,
    method = object$method
  )
  class(res) <- "summary.oneshot_frailty"
  res
}

#' @export
print.summary.oneshot_frailty <- function(x, digits = 4, signif.stars = TRUE, ...) {
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  cat("===================================================================\n")
  cat(sprintf(" One-Shot Device Testing Analysis with %s Baseline\n",
              if (x$baseline == "weibull") "Weibull" else "Exponential"))
  cat("===================================================================\n")
  cat(sprintf("Number of Components (M): %d (%s)\n", x$M, paste(x$components, collapse = ", ")))
  cat(sprintf("Estimation Algorithm:     %s\n", toupper(x$method)))
  cat(sprintf("Total Sample Size:        %d devices\n", x$sample_size))
  cat(sprintf("Log-Likelihood:           %.*f (df = %d)\n", digits, x$loglik, x$n_params))
  cat(sprintf("AIC: %.*f | AICc: %.*f | BIC: %.*f | HQIC: %.*f\n\n",
              digits, x$aic, digits, x$aicc, digits, x$bic, digits, x$hqic))
  
  cat("Model Coefficients & Asymptotic Inference:\n")
  stats::printCoefmat(x$coefficients, digits = digits, signif.stars = signif.stars, P.values = TRUE, has.Pvalue = TRUE)
  
  cat("\nAsymptotic Confidence Intervals:\n")
  print(round(x$confint, digits))
  
  if (x$frailty) {
    cat("\nFrailty Dependence Measures:\n")
    cat(sprintf("  Frailty Variance (beta):  %.*f\n", digits, x$coefficients["beta", "Estimate"]))
    cat(sprintf("  Kendall's rank tau:       %.*f (95%% CI: [%.*f, %.*f])\n",
                digits, x$kendall_tau$tau, digits, x$kendall_tau$confint["LCI"], digits, x$kendall_tau$confint["UCI"]))
  }
  cat("===================================================================\n")
  invisible(x)
}

#' @export
coef.oneshot_frailty <- function(object, ...) {
  object$coefficients
}

#' @export
vcov.oneshot_frailty <- function(object, ...) {
  object$vcov
}

#' @export
confint.oneshot_frailty <- function(object, parm, level = 0.95, ...) {
  cf <- object$coefficients
  pnames <- names(cf)
  if (missing(parm)) parm <- pnames
  else if (is.numeric(parm)) parm <- pnames[parm]
  
  se <- object$se
  z_val <- stats::qnorm(1 - (1 - level) / 2)
  
  ci <- cbind(
    LCI = cf[parm] - z_val * se[parm],
    UCI = cf[parm] + z_val * se[parm]
  )
  if ("beta" %in% parm) {
    ci["beta", "LCI"] <- max(0, ci["beta", "LCI"])
  }
  ci
}

#' @export
logLik.oneshot_frailty <- function(object, ...) {
  val <- object$loglik
  attr(val, "df") <- object$n_params
  attr(val, "nobs") <- object$sample_size
  class(val) <- "logLik"
  val
}

#' @export
AIC.oneshot_frailty <- function(object, ..., k = 2) {
  -2 * object$loglik + k * object$n_params
}

#' @export
BIC.oneshot_frailty <- function(object, ...) {
  -2 * object$loglik + log(object$sample_size) * object$n_params
}

#' @export
fitted.oneshot_frailty <- function(object, type = c("probs", "counts"), ...) {
  type <- match.arg(type)
  if (type == "probs") object$fitted_probs
  else object$fitted_counts
}

#' @export
residuals.oneshot_frailty <- function(object, type = c("pearson", "raw"), ...) {
  type <- match.arg(type)
  if (type == "pearson") object$residuals
  else object$count_mat - object$fitted_counts
}

#' @export
predict.oneshot_frailty <- function(object, newdata = NULL,
                                    type = c("probs", "reliability", "moments"),
                                    t = NULL, stress = NULL, k = 1, ...) {
  type <- match.arg(type)
  
  if (type == "reliability") {
    if (is.null(t)) stop("Please specify time points 't' for reliability prediction.")
    if (is.null(stress) && !is.null(newdata)) {
      stress <- newdata[[as.character(object$call$stress %||% "stress")]]
    }
    return(reliability_k_out_of_m(object, t = t, stress = stress, k = k))
  }
  
  if (type == "moments") {
    return(component_moments(object, stress = stress))
  }
  
  if (is.null(newdata)) {
    return(object$fitted_probs)
  }
  
  time_var <- as.character(object$call$time %||% "time")
  if (!(time_var %in% colnames(newdata))) {
    stop(sprintf("Time variable '%s' not found in newdata.", time_var))
  }
  tau_new <- as.numeric(newdata[[time_var]])
  
  stress_var <- if (ncol(object$X_cov) > 1) colnames(object$X_cov)[-1] else NULL
  if (!is.null(stress_var)) {
    cov_df <- newdata[, stress_var, drop = FALSE]
    X_new <- stats::model.matrix(~ ., data = cov_df)
  } else {
    X_new <- matrix(1, nrow = nrow(newdata), ncol = 1)
  }
  
  compute_probabilities(object$coefficients, X_new, tau_new, object$M,
                        object$baseline, object$frailty)$prob_matrix
}
