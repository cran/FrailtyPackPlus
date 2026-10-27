#' S3 Methods for Additive Shared Frailty Models
#'
#' Standard methods for printing, summarizing, extracting coefficients,
#' covariance matrices, log-likelihood, information criteria, residuals,
#' and predictions from fitted additive shared frailty model objects.
#'
#' @param x,object An object of class \code{"add_frailty"}, \code{"add_frailty_bayes"},
#'   or \code{"add_frailty_mle"}.
#' @param digits Number of decimal digits to print (default 4).
#' @param type For \code{residuals}, type of residual (\code{"cox_snell"});
#'   for \code{predict}, type of prediction (\code{"lp"} for linear predictor,
#'   \code{"survival"} for bivariate survival probabilities).
#' @param newdata Optional new data frame for predictions.
#' @param times Optional numeric vector of length 2 specifying \code{c(t1, t2)} for survival prediction.
#' @param k Penalty parameter for AIC (default 2).
#' @param ... Additional arguments passed to methods.
#'
#' @return Formatted console output, matrices, numeric vectors, or summaries.
#'
#' @importFrom stats AIC BIC coef logLik vcov residuals predict model.frame model.matrix dgamma rgamma rnorm pnorm dnorm rexp var runif uniroot density acf plogis qlogis ks.test
#' @importFrom graphics plot lines points legend abline polygon grid par
#' @importFrom grDevices rgb
#' @importFrom utils data head
#'
#' @name s3_methods
#' @export
print.add_frailty <- function(x, digits = 4, ...) {
  cat("\nAdditive Shared Frailty Regression Model\n")
  cat("========================================\n")
  cat(sprintf("Frailty distribution  : %s\n", x$model_info$frailty_type))
  cat(sprintf("Baseline distribution : %s\n", x$model_info$baseline_dist))
  cat(sprintf("Number of subjects    : %d\n", x$model_info$n_obs))
  cat(sprintf("Estimated parameters  : %d\n\n", x$model_info$n_pars))

  cat("Coefficients / Estimates (Point Estimates):\n")
  print(round(x$coefficients, digits))
  cat("\n")
  cat(sprintf("Kendall's tau (tau)   : %.4f\n", x$kendall_tau))
  cat(sprintf("Frailty Variance      : %.4f\n", x$frailty_variance))
  cat(sprintf("Log-Likelihood        : %.4f\n", x$criteria$LogLik))
  cat(sprintf("AIC: %.2f | BIC: %.2f | AICc: %.2f | HQIC: %.2f",
              x$criteria$AIC, x$criteria$BIC, x$criteria$AICc, x$criteria$HQIC))
  if (!is.null(x$criteria$DIC) && !is.na(x$criteria$DIC)) {
    cat(sprintf(" | DIC: %.2f", x$criteria$DIC))
  }
  cat("\n\n")
  invisible(x)
}

#' @rdname s3_methods
#' @export
summary.add_frailty <- function(object, ...) {
  is_bayes <- inherits(object, "add_frailty_bayes")

  res <- list(
    call = object$model_info$call,
    frailty_type = object$model_info$frailty_type,
    baseline_dist = object$model_info$baseline_dist,
    is_bayes = is_bayes,
    n_obs = object$model_info$n_obs,
    n_pars = object$model_info$n_pars,
    estimates = object$estimates,
    diagnostics = if (is_bayes) object$diagnostics else NULL,
    criteria = object$criteria,
    kendall_tau = object$kendall_tau,
    frailty_variance = object$frailty_variance
  )
  class(res) <- "summary.add_frailty"
  res
}

#' @rdname s3_methods
#' @export
print.summary.add_frailty <- function(x, digits = 4, ...) {
  cat("\n=======================================================================\n")
  cat("Additive Shared Frailty Model Summary\n")
  cat("=======================================================================\n")
  cat(sprintf("Frailty Distribution  : %s\n", x$frailty_type))
  cat(sprintf("Baseline Distribution : %s\n", x$baseline_dist))
  cat(sprintf("Estimation Framework  : %s\n", if (x$is_bayes) "Bayesian MCMC" else "Maximum Likelihood"))
  cat(sprintf("Number of Observations: %d\n", x$n_obs))
  cat("-----------------------------------------------------------------------\n")

  cat("\nParameter Estimates:\n")
  est_print <- x$estimates
  num_cols <- vapply(est_print, is.numeric, logical(1))
  est_print[num_cols] <- lapply(est_print[num_cols], round, digits = digits)
  print(est_print, row.names = FALSE)

  if (!is.null(x$diagnostics)) {
    cat("\nMCMC Convergence Diagnostics:\n")
    diag_print <- x$diagnostics
    num_cols_d <- vapply(diag_print, is.numeric, logical(1))
    diag_print[num_cols_d] <- lapply(diag_print[num_cols_d], round, digits = digits)
    print(diag_print, row.names = FALSE)
  }

  cat("\nDependence and Heterogeneity:\n")
  cat(sprintf("  Kendall's tau (tau)   : %.4f\n", x$kendall_tau))
  cat(sprintf("  Theoretical Frailty Var: %.4f\n", x$frailty_variance))

  cat("\nModel Selection Criteria:\n")
  cat(sprintf("  Log-Likelihood : %.4f\n", x$criteria$LogLik))
  cat(sprintf("  AIC            : %.4f\n", x$criteria$AIC))
  cat(sprintf("  BIC            : %.4f\n", x$criteria$BIC))
  cat(sprintf("  AICc           : %.4f\n", x$criteria$AICc))
  cat(sprintf("  HQIC           : %.4f\n", x$criteria$HQIC))
  if (!is.null(x$criteria$DIC) && !is.na(x$criteria$DIC)) {
    cat(sprintf("  DIC            : %.4f\n", x$criteria$DIC))
  }
  cat("=======================================================================\n\n")
  invisible(x)
}

#' @rdname s3_methods
#' @export
coef.add_frailty <- function(object, ...) {
  object$coefficients
}

#' @rdname s3_methods
#' @export
vcov.add_frailty <- function(object, ...) {
  if (inherits(object, "add_frailty_bayes")) {
    stats::var(object$posterior_samples)
  } else {
    object$vcov
  }
}

#' @rdname s3_methods
#' @export
logLik.add_frailty <- function(object, ...) {
  ll <- object$criteria$LogLik
  attr(ll, "df") <- object$model_info$n_pars
  attr(ll, "nobs") <- object$model_info$n_obs
  class(ll) <- "logLik"
  ll
}

#' @rdname s3_methods
#' @export
AIC.add_frailty <- function(object, ..., k = 2) {
  object$criteria$AIC
}

#' @rdname s3_methods
#' @export
BIC.add_frailty <- function(object, ...) {
  object$criteria$BIC
}

#' @rdname s3_methods
#' @export
residuals.add_frailty <- function(object, type = c("cox_snell"), ...) {
  type <- match.arg(type)
  pars <- object$coefficients
  if (object$model_info$baseline_dist == "gw") {
    cumh_fn <- cumh_gw
  } else {
    cumh_fn <- cumh_gll2
  }

  zeta1 <- pars[1]; delta1 <- pars[2]; xi1 <- pars[3]
  zeta2 <- pars[4]; delta2 <- pars[5]; xi2 <- pars[6]

  data.frame(
    r1 = cumh_fn(1, delta1, zeta1, xi1),
    r2 = cumh_fn(1, delta2, zeta2, xi2)
  )
}

#' @rdname s3_methods
#' @export
predict.add_frailty <- function(object, newdata = NULL, type = c("lp", "survival"), times = NULL, ...) {
  type <- match.arg(type)
  pars <- object$coefficients
  n_pos <- if (object$model_info$frailty_type == "genlindley") 8 else 7
  n_pars <- length(pars)

  beta <- if (n_pars > n_pos) pars[(n_pos + 1):n_pars] else numeric(0)

  if (type == "lp") {
    if (length(beta) == 0) return(rep(0, if (is.null(newdata)) object$model_info$n_obs else nrow(newdata)))
    if (is.null(newdata)) {
      stop("'newdata' required for linear predictor calculation.")
    }
    X <- as.matrix(newdata)
    return(as.vector(X %*% beta))
  } else if (type == "survival") {
    if (is.null(times) || length(times) != 2) {
      stop("'times' must be a numeric vector of length 2: c(t1, t2).")
    }
    t1 <- times[1]
    t2 <- times[2]
    rho <- 1.0
    if (length(beta) > 0 && !is.null(newdata)) {
      X <- as.matrix(newdata)
      rho <- exp(mean(X %*% beta))
    }
    s_arg <- (t1 + t2) * rho
    if (object$model_info$baseline_dist == "gw") {
      Phi1 <- cumh_gw(t1, pars[2], pars[1], pars[3])
      Phi2 <- cumh_gw(t2, pars[5], pars[4], pars[6])
    } else {
      Phi1 <- cumh_gll2(t1, pars[2], pars[1], pars[3])
      Phi2 <- cumh_gll2(t2, pars[5], pars[4], pars[6])
    }

    if (object$model_info$frailty_type == "genlindley") {
      L0 <- laplace_genlindley(s_arg, pars[7], pars[8])
    } else if (object$model_info$frailty_type == "wlindley") {
      L0 <- laplace_wlindley(s_arg, pars[7])
    } else if (object$model_info$frailty_type == "gamma") {
      L0 <- laplace_gamma(s_arg, pars[7])
    } else {
      L0 <- laplace_invgaussian(s_arg, pars[7])
    }
    return(exp(-(Phi1 + Phi2)) * L0)
  }
}
