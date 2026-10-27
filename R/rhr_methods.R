# S3 Generic Methods for Shared RHR Frailty Model Fits

#' @importFrom stats AIC BIC coef logLik vcov
#' @export
print.rhr_bayes_fit <- function(x, ...) {
  cat("\n================================================================\n")
  cat(" Bayesian Shared RHR Frailty Model Fit\n")
  cat("================================================================\n")
  cat("Frailty Distribution:  ", x$model_info$frailty, "\n")
  cat("Baseline Distribution: ", x$model_info$baseline, "\n")
  cat("Number of Pairs (n):   ", x$model_info$n_obs, "\n")
  cat("MCMC Iterations:       ", x$model_info$n_iter, " (Burn-in: ", x$model_info$burn_in, ", Thin: ", x$model_info$thin, ")\n")
  cat("Kendall's Tau (tau):   ", round(x$kendall_tau, 4), "\n")
  cat("Frailty Variance V(W): ", round(x$frailty_variance, 4), "\n")
  cat("\nPosterior Estimates (SELF / Posterior Mean):\n")
  print(round(x$estimates$SELF, 5))
  cat("\nInformation Criteria:\n")
  cat("  AIC:", round(x$criteria$AIC, 2), "  BIC:", round(x$criteria$BIC, 2),
      "  DIC:", round(x$criteria$DIC, 2), "  LPML:", round(x$criteria$LPML, 2), "\n")
  cat("================================================================\n")
  invisible(x)
}

#' @export
summary.rhr_bayes_fit <- function(object, ...) {
  out <- list(
    summary_table = object$summary,
    model_info = object$model_info,
    criteria = object$criteria,
    kendall_tau = object$kendall_tau,
    frailty_variance = object$frailty_variance
  )
  class(out) <- "summary.rhr_bayes_fit"
  return(out)
}

#' @export
print.summary.rhr_bayes_fit <- function(x, ...) {
  cat("\n=======================================================================\n")
  cat(" Summary: Bayesian Shared RHR Frailty Model\n")
  cat("=======================================================================\n")
  cat("Frailty: ", x$model_info$frailty, " | Baseline: ", x$model_info$baseline,
      " | Sample Size: ", x$model_info$n_obs, "\n\n")
  cat("Posterior Parameter Estimates under Multiple Loss Functions:\n")
  print(x$summary_table, row.names = FALSE)
  cat("\nDependence and Heterogeneity:\n")
  cat("  Kendall's Tau:      ", round(x$kendall_tau, 5), "\n")
  cat("  Frailty Variance:   ", round(x$frailty_variance, 5), "\n")
  cat("\nModel Comparison Metrics:\n")
  cat("  AIC: ", round(x$criteria$AIC, 2), " | BIC: ", round(x$criteria$BIC, 2),
      " | AICc: ", round(x$criteria$AICc, 2), " | HQIC: ", round(x$criteria$HQIC, 2), "\n")
  cat("  DIC: ", round(x$criteria$DIC, 2), " | pD:  ", round(x$criteria$pD, 2),
      " | EAIC: ", round(x$criteria$EAIC, 2), " | EBIC: ", round(x$criteria$EBIC, 2), "\n")
  cat("  LPML (CPO): ", round(x$criteria$LPML, 2), "\n")
  cat("=======================================================================\n")
  invisible(x)
}

#' @export
coef.rhr_bayes_fit <- function(object, loss = c("SELF", "ELF", "PLF"), ...) {
  loss <- match.arg(loss)
  object$estimates[[loss]]
}

#' @export
logLik.rhr_bayes_fit <- function(object, ...) {
  val <- -0.5 * object$criteria$dev_at_mean
  attr(val, "df") <- length(object$estimates$SELF)
  attr(val, "nobs") <- object$model_info$n_obs
  class(val) <- "logLik"
  return(val)
}

#' @export
AIC.rhr_bayes_fit <- function(object, ...) {
  object$criteria$AIC
}

#' @export
BIC.rhr_bayes_fit <- function(object, ...) {
  object$criteria$BIC
}

#' @export
print.rhr_mle_fit <- function(x, ...) {
  cat("\n================================================================\n")
  cat(" Maximum Likelihood Shared RHR Frailty Model Fit\n")
  cat("================================================================\n")
  cat("Frailty Distribution:  ", x$model_info$frailty, "\n")
  cat("Baseline Distribution: ", x$model_info$baseline, "\n")
  cat("Number of Pairs (n):   ", x$model_info$n_obs, "\n")
  cat("Log-Likelihood:        ", round(x$loglik, 4), "\n")
  cat("Kendall's Tau (tau):   ", round(x$kendall_tau, 4), "\n")
  cat("Frailty Variance V(W): ", round(x$frailty_variance, 4), "\n")
  cat("\nCoefficients:\n")
  print(round(x$coefficients, 5))
  cat("\nInformation Criteria:\n")
  cat("  AIC:", round(x$criteria$AIC, 2), "  BIC:", round(x$criteria$BIC, 2),
      "  AICc:", round(x$criteria$AICc, 2), "  HQIC:", round(x$criteria$HQIC, 2), "\n")
  cat("================================================================\n")
  invisible(x)
}

#' @export
summary.rhr_mle_fit <- function(object, ...) {
  out <- list(
    summary_table = object$summary,
    model_info = object$model_info,
    criteria = object$criteria,
    loglik = object$loglik,
    kendall_tau = object$kendall_tau,
    frailty_variance = object$frailty_variance
  )
  class(out) <- "summary.rhr_mle_fit"
  return(out)
}

#' @export
print.summary.rhr_mle_fit <- function(x, ...) {
  cat("\n=======================================================================\n")
  cat(" Summary: Maximum Likelihood Shared RHR Frailty Model\n")
  cat("=======================================================================\n")
  cat("Frailty: ", x$model_info$frailty, " | Baseline: ", x$model_info$baseline,
      " | Sample Size: ", x$model_info$n_obs, "\n\n")
  cat("Parameter Estimates and Wald Confidence Intervals:\n")
  print(x$summary_table, row.names = FALSE)
  cat("\nDependence and Heterogeneity:\n")
  cat("  Kendall's Tau:      ", round(x$kendall_tau, 5), "\n")
  cat("  Frailty Variance:   ", round(x$frailty_variance, 5), "\n")
  cat("\nModel Fit and Information Criteria:\n")
  cat("  Log-Likelihood: ", round(x$loglik, 2), "\n")
  cat("  AIC: ", round(x$criteria$AIC, 2), " | BIC: ", round(x$criteria$BIC, 2),
      " | AICc: ", round(x$criteria$AICc, 2), " | HQIC: ", round(x$criteria$HQIC, 2), "\n")
  cat("=======================================================================\n")
  invisible(x)
}

#' @export
coef.rhr_mle_fit <- function(object, ...) {
  object$coefficients
}

#' @export
vcov.rhr_mle_fit <- function(object, ...) {
  object$vcov
}

#' @export
logLik.rhr_mle_fit <- function(object, ...) {
  val <- object$loglik
  attr(val, "df") <- length(object$coefficients)
  attr(val, "nobs") <- object$model_info$n_obs
  class(val) <- "logLik"
  return(val)
}

#' @export
AIC.rhr_mle_fit <- function(object, ...) {
  object$criteria$AIC
}

#' @export
BIC.rhr_mle_fit <- function(object, ...) {
  object$criteria$BIC
}
