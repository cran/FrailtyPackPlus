#' Methods for ziprom_fit Objects
#'
#' Standard S3 methods for inspecting, summarizing, and extracting components
#' from fitted zero-inflated promotion cure-rate models.
#'
#' @param x An object of class \code{"ziprom_fit"} or \code{"summary.ziprom_fit"}.
#' @param object An object of class \code{"ziprom_fit"}.
#' @param digits Integer indicating the number of decimal digits to print. Default is 4.
#' @param level Confidence level for intervals (e.g. 0.90, 0.95, 0.99). Default is 0.95.
#' @param parm Optional vector of parameter names or indices for \code{confint}.
#' @param k Penalty parameter for \code{AIC} (default \code{k = 2}).
#' @param type Character string indicating residual type: \code{"quantile"} (default),
#'   \code{"cox-snell"}, \code{"deviance"}, or \code{"response"}.
#' @param ... Additional arguments passed to or from other methods.
#' @importFrom stats AIC BIC coef confint logLik residuals vcov
#'
#' @return
#' \itemize{
#'   \item \code{print.ziprom_fit}: Invisible copy of the model object.
#'   \item \code{summary.ziprom_fit}: A summary object of class \code{"summary.ziprom_fit"}.
#'   \item \code{coef.ziprom_fit}: Named vector of estimated model coefficients.
#'   \item \code{vcov.ziprom_fit}: Asymptotic variance-covariance matrix.
#'   \item \code{confint.ziprom_fit}: Matrix of confidence intervals.
#'   \item \code{logLik.ziprom_fit}: Log-likelihood value with attributes \code{df} and \code{nobs}.
#'   \item \code{AIC.ziprom_fit}: Akaike Information Criterion.
#'   \item \code{BIC.ziprom_fit}: Bayesian Information Criterion.
#'   \item \code{residuals.ziprom_fit}: Vector of residuals.
#' }
#'
#' @export
print.ziprom_fit <- function(x, digits = 4, ...) {
  cat("\n=======================================================\n")
  cat("Zero-Inflated Promotion Cure-Rate Model with Frailty\n")
  cat("=======================================================\n")
  cat("Activation Scheme :", x$activation, "\n")
  cat("Baseline Hazard   :", x$baseline, "\n")
  cat("Gamma Frailty     :", ifelse(x$frailty, "Yes", "No"), "\n")
  cat("Sample Size (n)   :", x$n, "\n")
  cat("Zero Failures     :", x$n_zeros, sprintf("(%.1f%%)", 100 * x$n_zeros / x$n), "\n")
  cat("Observed Events   :", x$n_events, sprintf("(%.1f%%)", 100 * x$n_events / x$n), "\n")
  cat("Right-Censored    :", x$n_censored, sprintf("(%.1f%%)", 100 * x$n_censored / x$n), "\n")
  cat("-------------------------------------------------------\n")
  cat("Log-Likelihood    :", round(x$loglik, digits), "\n")
  cat("AIC               :", round(x$aic, digits), "\n")
  cat("BIC               :", round(x$bic, digits), "\n")
  cat("-------------------------------------------------------\n")
  cat("Estimated Coefficients:\n")
  print(round(x$coefficients, digits))
  if (x$frailty) {
    cat(sprintf("\nEstimated Frailty Variance (alpha): %.4f (SE: %.4f)\n",
                x$coef_list$alpha, x$coef_list$alpha_se))
  }
  cat(sprintf("Estimated Mean p0 (zeros): %.4f (SE: %.4f)\n", x$p0_mean, x$p0_mean_se))
  cat(sprintf("Estimated Mean p1 (cure) : %.4f (SE: %.4f)\n", x$p1_mean, x$p1_mean_se))
  cat("=======================================================\n\n")
  invisible(x)
}

#' @rdname print.ziprom_fit
#' @export
summary.ziprom_fit <- function(object, level = 0.95, ...) {
  z_crit <- stats::qnorm(1 - (1 - level) / 2)

  ci_mat <- cbind(
    Lower = object$coefficients - z_crit * object$se,
    Upper = object$coefficients + z_crit * object$se
  )

  tab <- data.frame(
    Estimate = object$coefficients,
    `Std. Error` = object$se,
    `z value` = object$z_values,
    `Pr(>|z|)` = object$p_values,
    `CI Lower` = ci_mat[, 1],
    `CI Upper` = ci_mat[, 2],
    check.names = FALSE
  )

  res <- list(
    call = object$call,
    activation = object$activation,
    baseline = object$baseline,
    frailty = object$frailty,
    n = object$n,
    n_zeros = object$n_zeros,
    n_events = object$n_events,
    n_censored = object$n_censored,
    loglik = object$loglik,
    aic = object$aic,
    bic = object$bic,
    hqic = object$hqic,
    caic = object$caic,
    coef_table = tab,
    coef_list = object$coef_list,
    p0_mean = object$p0_mean,
    p0_mean_se = object$p0_mean_se,
    p1_mean = object$p1_mean,
    p1_mean_se = object$p1_mean_se,
    level = level
  )
  class(res) <- "summary.ziprom_fit"
  res
}

#' @rdname print.ziprom_fit
#' @export
print.summary.ziprom_fit <- function(x, digits = 4, ...) {
  cat("\n========================================================================\n")
  cat("Summary of Zero-Inflated Promotion Cure-Rate Model with Frailty\n")
  cat("========================================================================\n")
  cat("Activation Scheme :", x$activation, "\n")
  cat("Baseline Hazard   :", x$baseline, "\n")
  cat("Gamma Frailty     :", ifelse(x$frailty, "Yes", "No"), "\n")
  cat("Sample Size (n)   :", x$n, "\n")
  cat("Zeros: ", x$n_zeros, " | Events: ", x$n_events, " | Censored: ", x$n_censored, "\n", sep = "")
  cat("------------------------------------------------------------------------\n")
  cat(sprintf("Log-Lik: %.4f | AIC: %.4f | BIC: %.4f | HQIC: %.4f | CAIC: %.4f\n",
              x$loglik, x$aic, x$bic, x$hqic, x$caic))
  cat("------------------------------------------------------------------------\n")
  cat(sprintf("Coefficients Table (%d%% Confidence Intervals):\n\n", round(x$level * 100)))

  formatted_tab <- data.frame(
    Estimate = sprintf("%.*f", digits, x$coef_table$Estimate),
    `Std. Error` = sprintf("%.*f", digits, x$coef_table$`Std. Error`),
    `z value` = sprintf("%.*f", digits, x$coef_table$`z value`),
    `Pr(>|z|)` = format.pval(x$coef_table$`Pr(>|z|)`, digits = digits, eps = 1e-4),
    `CI Lower` = sprintf("%.*f", digits, x$coef_table$`CI Lower`),
    `CI Upper` = sprintf("%.*f", digits, x$coef_table$`CI Upper`),
    row.names = rownames(x$coef_table),
    check.names = FALSE
  )
  print(formatted_tab, quote = FALSE, right = TRUE)

  cat("------------------------------------------------------------------------\n")
  if (x$frailty) {
    cat(sprintf("Frailty Variance (alpha): %.4f (SE: %.4f)\n",
                x$coef_list$alpha, x$coef_list$alpha_se))
  }
  cat(sprintf("Zero Proportion (p0)    : %.4f (SE: %.4f)\n", x$p0_mean, x$p0_mean_se))
  cat(sprintf("Cure Proportion (p1)    : %.4f (SE: %.4f)\n", x$p1_mean, x$p1_mean_se))
  cat("========================================================================\n\n")
  invisible(x)
}

#' @rdname print.ziprom_fit
#' @export
coef.ziprom_fit <- function(object, ...) {
  object$coefficients
}

#' @rdname print.ziprom_fit
#' @export
vcov.ziprom_fit <- function(object, ...) {
  object$vcov
}

#' @rdname print.ziprom_fit
#' @export
confint.ziprom_fit <- function(object, parm, level = 0.95, ...) {
  cf <- object$coefficients
  pnames <- names(cf)

  if (missing(parm)) {
    parm <- pnames
  } else if (is.numeric(parm)) {
    parm <- pnames[parm]
  }

  z_crit <- stats::qnorm(1 - (1 - level) / 2)
  se <- object$se[parm]
  est <- cf[parm]

  pct <- paste(round(100 * c((1 - level) / 2, 1 - (1 - level) / 2), 1), "%")
  ci <- cbind(est - z_crit * se, est + z_crit * se)
  colnames(ci) <- pct
  rownames(ci) <- parm
  ci
}

#' @rdname print.ziprom_fit
#' @export
logLik.ziprom_fit <- function(object, ...) {
  structure(
    object$loglik,
    df = object$npar,
    nobs = object$n,
    class = "logLik"
  )
}

#' @rdname print.ziprom_fit
#' @export
AIC.ziprom_fit <- function(object, ..., k = 2) {
  -2 * object$loglik + k * object$npar
}

#' @rdname print.ziprom_fit
#' @export
BIC.ziprom_fit <- function(object, ...) {
  -2 * object$loglik + log(object$n) * object$npar
}

#' @rdname print.ziprom_fit
#' @export
residuals.ziprom_fit <- function(object, type = c("quantile", "cox-snell", "deviance", "response"), ...) {
  type <- match.arg(type)

  if (type == "quantile") {
    return(object$residuals)
  }

  time_vec <- object$time
  status_vec <- object$status
  alpha_val <- if (object$frailty) object$coef_list$alpha else 0

  shape_vec <- exp(as.vector(object$X3 %*% object$coef_list$beta3))
  scale_vec <- exp(as.vector(object$X4 %*% object$coef_list$beta4))
  probs <- compute_p0_p1(object$X1, object$X2, object$coef_list$beta1, object$coef_list$beta2)

  s_star <- ifelse(time_vec == 0, 1,
                   surv_proper(time_vec, probs$theta, alpha_val, shape_vec, scale_vec,
                               object$activation, object$baseline))
  s_pop <- ifelse(time_vec == 0, 1 - probs$p0, probs$p1 + probs$susceptible * s_star)
  s_pop <- pmin(pmax(s_pop, 1e-12), 1 - 1e-12)

  if (type == "cox-snell") {
    -log(s_pop)
  } else if (type == "deviance") {
    cs <- -log(s_pop)
    sign(status_vec - cs) * sqrt(pmax(2 * (cs - status_vec - status_vec * log(pmax(cs, 1e-12))), 0))
  } else if (type == "response") {
    status_vec - (1 - s_pop)
  }
}
