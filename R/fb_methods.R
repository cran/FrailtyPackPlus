#' Print Methods for Frailty Models
#'
#' @return The evaluated object \code{x}, invisibly.
#' @param x An object of class \code{"frailty_bayes"} or \code{"frailty_mle"}.
#' @param ... Additional arguments.
#' @export
print.frailty_bayes <- function(x, ...) {
  cat("\n=======================================================\n")
  cat("  Bayesian Frailty Survival Regression Model\n")
  cat("=======================================================\n")
  cat("Baseline Distribution :", x$baseline, "\n")
  cat("Frailty Distribution  :", x$frailty, "\n")
  cat("Model Formulation     :", x$type, "\n")
  cat("Number of Observations:", x$data$n, "\n")
  cat("MCMC Iterations       :", nrow(x$chain), "(post burn-in)\n")
  cat("Acceptance Rate       :", round(x$accept_rate * 100, 2), "%\n")
  cat("-------------------------------------------------------\n")
  cat("Model Fit Information:\n")
  cat("  Log-Likelihood :", round(x$diagnostics$LogLik, 3), "\n")
  cat("  DIC            :", round(x$diagnostics$DIC, 3), "(p_D =", round(x$diagnostics$p_D, 3), ")\n")
  cat("  EAIC           :", round(x$diagnostics$EAIC, 3), "\n")
  cat("  EBIC           :", round(x$diagnostics$EBIC, 3), "\n")
  cat("  LPML           :", round(x$diagnostics$LPML, 3), "\n")
  cat("=======================================================\n\n")
  invisible(x)
}

#' @export
print.frailty_mle <- function(x, ...) {
  cat("\n=======================================================\n")
  cat("  Maximum Likelihood Frailty Regression Model\n")
  cat("=======================================================\n")
  cat("Baseline Distribution :", x$baseline, "\n")
  cat("Frailty Distribution  :", x$frailty, "\n")
  cat("Model Formulation     :", x$type, "\n")
  cat("Number of Observations:", x$n_obs, "\n")
  cat("Log-Likelihood        :", round(x$logLik, 3), "\n")
  cat("AIC                   :", round(x$AIC, 3), "\n")
  cat("BIC                   :", round(x$BIC, 3), "\n")
  cat("=======================================================\n\n")
  invisible(x)
}

#' Summary Method for Bayesian Frailty Models
#'
#' @return An object of class \code{"summary.frailty_bayes"}.
#' @param object Object of class \code{"frailty_bayes"}.
#' @param ... Additional arguments.
#' @export
summary.frailty_bayes <- function(object, ...) {
  structure(
    list(
      call = object$call,
      baseline = object$baseline,
      frailty = object$frailty,
      type = object$type,
      n_obs = object$data$n,
      n_iter = nrow(object$chain),
      accept_rate = object$accept_rate,
      coefficients = object$summary,
      diagnostics = object$diagnostics
    ),
    class = "summary.frailty_bayes"
  )
}

#' @export
print.summary.frailty_bayes <- function(x, ...) {
  cat("\n===============================================================================\n")
  cat("  Bayesian Frailty Model: Posterior Estimates & Diagnostics\n")
  cat("===============================================================================\n")
  cat("Baseline:", x$baseline, "| Frailty:", x$frailty, "| Type:", x$type, "\n")
  cat("Observations:", x$n_obs, "| Saved Iterations:", x$n_iter, "| Acceptance:", round(x$accept_rate * 100, 2), "%\n")
  cat("-------------------------------------------------------------------------------\n")
  
  # Format table
  df <- x$coefficients
  num_cols <- c("Mean", "SD", "Median", "MCSE", "HPD_Lower", "HPD_Upper", "Quantile_2.5", "Quantile_97.5", "Geweke_Z")
  for (col in num_cols) {
    if (col %in% names(df)) df[[col]] <- round(df[[col]], 4)
  }
  print(df, row.names = FALSE)
  
  cat("-------------------------------------------------------------------------------\n")
  cat("Goodness-of-Fit / Comparison Criteria:\n")
  cat("  DIC  :", round(x$diagnostics$DIC, 3), "  | D_bar:", round(x$diagnostics$D_bar, 3),
      "  | D_hat:", round(x$diagnostics$D_hat, 3), "  | p_D:", round(x$diagnostics$p_D, 3), "\n")
  cat("  EAIC :", round(x$diagnostics$EAIC, 3), "  | EBIC :", round(x$diagnostics$EBIC, 3),
      "  | LPML :", round(x$diagnostics$LPML, 3), "\n")
  cat("===============================================================================\n\n")
  invisible(x)
}

#' Summary Method for MLE Frailty Models
#'
#' @return An object of class \code{"summary.frailty_mle"}.
#' @param object Object of class \code{"frailty_mle"}.
#' @param ... Additional arguments.
#' @export
summary.frailty_mle <- function(object, ...) {
  structure(
    list(
      call = object$call,
      baseline = object$baseline,
      frailty = object$frailty,
      type = object$type,
      n_obs = object$n_obs,
      coefficients = object$summary,
      logLik = object$logLik,
      AIC = object$AIC,
      BIC = object$BIC
    ),
    class = "summary.frailty_mle"
  )
}

#' @export
print.summary.frailty_mle <- function(x, ...) {
  cat("\n===============================================================================\n")
  cat("  Maximum Likelihood Frailty Model Summary\n")
  cat("===============================================================================\n")
  cat("Baseline:", x$baseline, "| Frailty:", x$frailty, "| Type:", x$type, "\n")
  cat("Observations:", x$n_obs, "| Log-Likelihood:", round(x$logLik, 3), "\n")
  cat("AIC:", round(x$AIC, 3), "| BIC:", round(x$BIC, 3), "\n")
  cat("-------------------------------------------------------------------------------\n")
  print(round(x$coefficients, 4))
  cat("===============================================================================\n\n")
  invisible(x)
}

#' S3 Plot Method for Bayesian Frailty Models
#'
#' @return Invisibly returns \code{NULL}.
#'
#' Produces trace plots, posterior density plots, autocorrelation plots,
#' estimated survival curves, or estimated hazard rate functions.
#'
#' @param x An object of class \code{"frailty_bayes"}.
#' @param which Character vector specifying plot types: \code{"trace"},
#'   \code{"density"}, \code{"autocorr"}, \code{"survival"}, or \code{"hazard"}.
#' @param ... Additional graphical parameters.
#' @export
plot.frailty_bayes <- function(x, which = c("trace", "density", "autocorr", "survival", "hazard"), ...) {
  which <- match.arg(which)
  chain <- x$chain
  p_names <- colnames(chain)
  n_p <- ncol(chain)
  
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par))
  
  if (which == "trace") {
    n_rows <- ceiling(sqrt(n_p))
    n_cols <- ceiling(n_p / n_rows)
    graphics::par(mfrow = c(n_rows, n_cols), mar = c(3, 3, 2, 1), mgp = c(1.8, 0.6, 0))
    for (j in 1:n_p) {
      graphics::plot(chain[, j], type = "l", col = "steelblue",
                     main = paste("Trace:", p_names[j]),
                     xlab = "Iteration", ylab = "Value")
      graphics::abline(h = mean(chain[, j]), col = "red", lty = 2)
    }
  } else if (which == "density") {
    n_rows <- ceiling(sqrt(n_p))
    n_cols <- ceiling(n_p / n_rows)
    graphics::par(mfrow = c(n_rows, n_cols), mar = c(3, 3, 2, 1), mgp = c(1.8, 0.6, 0))
    for (j in 1:n_p) {
      d <- stats::density(chain[, j])
      graphics::plot(d, col = "darkgreen", lwd = 2,
                     main = paste("Density:", p_names[j]),
                     xlab = "Value", ylab = "Posterior Density")
      graphics::polygon(d, col = grDevices::rgb(0, 0.5, 0, 0.2), border = NA)
      graphics::abline(v = mean(chain[, j]), col = "red", lty = 2)
    }
  } else if (which == "autocorr") {
    n_rows <- ceiling(sqrt(n_p))
    n_cols <- ceiling(n_p / n_rows)
    graphics::par(mfrow = c(n_rows, n_cols), mar = c(3, 3, 2, 1), mgp = c(1.8, 0.6, 0))
    for (j in 1:n_p) {
      stats::acf(chain[, j], main = paste("ACF:", p_names[j]), lag.max = 40)
    }
  } else if (which == "survival") {
    # Plot estimated population survival curve
    grid_t <- seq(0.1, 100, length.out = 200)
    if (grepl("bivariate", x$type)) {
      grid_t <- seq(0.1, max(c(x$data$time1, x$data$time2), na.rm = TRUE) * 1.1, length.out = 200)
    } else {
      grid_t <- seq(0.1, max(x$data$time, na.rm = TRUE) * 1.1, length.out = 200)
    }
    
    mean_p <- colMeans(chain)
    base_p <- mean_p[1:x$p_base]
    frail_p <- mean_p[(x$p_base + x$p_cov + 1):(x$p_base + x$p_cov + x$p_frail)]
    
    b_eval <- eval_baseline(grid_t, x$baseline, base_p)
    fr_eval <- eval_laplace(b_eval$H0, x$frailty, frail_p)
    
    graphics::plot(grid_t, fr_eval$L0, type = "l", col = "blue", lwd = 2.5,
                   ylim = c(0, 1), xlab = "Time", ylab = "Unconditional Survival S(t)",
                   main = paste("Estimated Marginal Survival Curve (", x$baseline, "+", x$frailty, ")"))
    graphics::lines(grid_t, b_eval$S0, col = "darkgray", lty = 2, lwd = 1.8)
    graphics::legend("topright", legend = c("Marginal (Frailty)", "Baseline S0(t)"),
                     col = c("blue", "darkgray"), lty = c(1, 2), lwd = c(2.5, 1.8))
  } else if (which == "hazard") {
    grid_t <- seq(0.1, 100, length.out = 200)
    if (grepl("bivariate", x$type)) {
      grid_t <- seq(0.1, max(c(x$data$time1, x$data$time2), na.rm = TRUE) * 1.1, length.out = 200)
    } else {
      grid_t <- seq(0.1, max(x$data$time, na.rm = TRUE) * 1.1, length.out = 200)
    }
    mean_p <- colMeans(chain)
    base_p <- mean_p[1:x$p_base]
    b_eval <- eval_baseline(grid_t, x$baseline, base_p)
    graphics::plot(grid_t, b_eval$h0, type = "l", col = "firebrick", lwd = 2.5,
                   xlab = "Time", ylab = "Baseline Hazard Rate h0(t)",
                   main = paste("Estimated Baseline Hazard Curve (", x$baseline, ")"))
  }
}

#' Model Extraction Methods
#'
#' @return The extracted model quantity.
#' @param object Fitted frailty model object.
#' @param ... Additional arguments.
#' @export
logLik.frailty_bayes <- function(object, ...) {
  val <- object$diagnostics$LogLik
  attr(val, "df") <- object$diagnostics$n_params
  attr(val, "nobs") <- object$diagnostics$n_obs
  class(val) <- "logLik"
  val
}

#' @export
logLik.frailty_mle <- function(object, ...) {
  val <- object$logLik
  attr(val, "df") <- object$n_params
  attr(val, "nobs") <- object$n_obs
  class(val) <- "logLik"
  val
}

#' @export
AIC.frailty_bayes <- function(object, ..., k = 2) {
  object$diagnostics$EAIC
}

#' @export
AIC.frailty_mle <- function(object, ..., k = 2) {
  object$AIC
}

#' @export
BIC.frailty_bayes <- function(object, ...) {
  object$diagnostics$EBIC
}

#' @export
BIC.frailty_mle <- function(object, ...) {
  object$BIC
}

#' @export
coef.frailty_bayes <- function(object, ...) {
  colMeans(object$chain)
}

#' @export
coef.frailty_mle <- function(object, ...) {
  object$coefficients
}

#' @export
vcov.frailty_bayes <- function(object, ...) {
  stats::cov(object$chain)
}

#' @export
vcov.frailty_mle <- function(object, ...) {
  object$vcov
}

#' Residuals for Frailty Models
#'
#' Computes Cox-Snell, Martingale, or deviance residuals.
#'
#' @param object Fitted frailty model object.
#' @param type Type of residual: \code{"cox-snell"}, \code{"martingale"}, or \code{"deviance"}.
#' @param ... Additional arguments.
#' @return A numeric vector or matrix of residuals.
#' @export
residuals.frailty_bayes <- function(object, type = c("cox-snell", "martingale", "deviance"), ...) {
  type <- match.arg(type)
  mean_p <- colMeans(object$chain)
  base_p <- mean_p[1:object$p_base]
  beta <- if (object$p_cov > 0) mean_p[(object$p_base + 1):(object$p_base + object$p_cov)] else numeric(0)
  
  if (grepl("bivariate", object$type)) {
    b1 <- eval_baseline(object$data$time1, object$baseline, base_p)
    b2 <- eval_baseline(object$data$time2, object$baseline, base_p)
    lp1 <- if (object$p_cov > 0) as.vector(object$data$X1 %*% beta) else rep(0, object$data$n)
    lp2 <- if (object$p_cov > 0) as.vector(object$data$X2 %*% beta) else rep(0, object$data$n)
    
    H1 <- b1$H0 * exp(lp1)
    H2 <- b2$H0 * exp(lp2)
    
    if (type == "cox-snell") {
      return(data.frame(res1 = H1, res2 = H2))
    } else if (type == "martingale") {
      return(data.frame(res1 = object$data$status1 - H1, res2 = object$data$status2 - H2))
    } else {
      m1 <- object$data$status1 - H1
      m2 <- object$data$status2 - H2
      d1 <- sign(m1) * sqrt(pmax(0, -2 * (m1 + ifelse(object$data$status1 == 1, log(pmax(H1, 1e-12)), 0))))
      d2 <- sign(m2) * sqrt(pmax(0, -2 * (m2 + ifelse(object$data$status2 == 1, log(pmax(H2, 1e-12)), 0))))
      return(data.frame(res1 = d1, res2 = d2))
    }
  } else {
    b <- eval_baseline(object$data$time, object$baseline, base_p)
    lp <- if (object$p_cov > 0) as.vector(object$data$X %*% beta) else rep(0, object$data$n)
    H <- b$H0 * exp(lp)
    
    if (type == "cox-snell") {
      return(H)
    } else if (type == "martingale") {
      return(object$data$status - H)
    } else {
      m <- object$data$status - H
      return(sign(m) * sqrt(pmax(0, -2 * (m + ifelse(object$data$status == 1, log(pmax(H, 1e-12)), 0)))))
    }
  }
}

#' @export
residuals.frailty_mle <- function(object, type = c("cox-snell", "martingale", "deviance"), ...) {
  type <- match.arg(type)
  est_p <- object$coefficients
  n_base <- get_n_base_params(object$baseline)
  base_p <- est_p[1:n_base]
  beta <- if (object$data$p_cov > 0) est_p[(n_base + 1):(n_base + object$data$p_cov)] else numeric(0)
  
  if (grepl("bivariate", object$type)) {
    b1 <- eval_baseline(object$data$time1, object$baseline, base_p)
    b2 <- eval_baseline(object$data$time2, object$baseline, base_p)
    lp1 <- if (object$data$p_cov > 0) as.vector(object$data$X1 %*% beta) else rep(0, object$data$n)
    lp2 <- if (object$data$p_cov > 0) as.vector(object$data$X2 %*% beta) else rep(0, object$data$n)
    H1 <- b1$H0 * exp(lp1)
    H2 <- b2$H0 * exp(lp2)
    if (type == "cox-snell") return(data.frame(res1 = H1, res2 = H2))
    if (type == "martingale") return(data.frame(res1 = object$data$status1 - H1, res2 = object$data$status2 - H2))
    m1 <- object$data$status1 - H1
    m2 <- object$data$status2 - H2
    return(data.frame(
      res1 = sign(m1) * sqrt(pmax(0, -2 * (m1 + ifelse(object$data$status1 == 1, log(pmax(H1, 1e-12)), 0)))),
      res2 = sign(m2) * sqrt(pmax(0, -2 * (m2 + ifelse(object$data$status2 == 1, log(pmax(H2, 1e-12)), 0))))
    ))
  } else {
    b <- eval_baseline(object$data$time, object$baseline, base_p)
    lp <- if (object$data$p_cov > 0) as.vector(object$data$X %*% beta) else rep(0, object$data$n)
    H <- b$H0 * exp(lp)
    if (type == "cox-snell") return(H)
    if (type == "martingale") return(object$data$status - H)
    m <- object$data$status - H
    return(sign(m) * sqrt(pmax(0, -2 * (m + ifelse(object$data$status == 1, log(pmax(H, 1e-12)), 0)))))
  }
}

#' Predict Method for Frailty Models
#'
#' Predicts survival probabilities, cumulative hazards, or frailty estimates.
#'
#' @param object Fitted model object.
#' @param newdata New data frame.
#' @param times Time points for evaluation.
#' @param type Prediction type: \code{"survival"}, \code{"cumhaz"}, or \code{"hazard"}.
#' @param ... Additional arguments.
#' @return A numeric vector or matrix of predictions.
#' @export
predict.frailty_bayes <- function(object, newdata = NULL, times = NULL,
                                 type = c("survival", "cumhaz", "hazard"), ...) {
  type <- match.arg(type)
  if (is.null(times)) times <- seq(1, 100, length.out = 50)
  
  mean_p <- colMeans(object$chain)
  base_p <- mean_p[1:object$p_base]
  frail_p <- mean_p[(object$p_base + object$p_cov + 1):(object$p_base + object$p_cov + object$p_frail)]
  
  b_eval <- eval_baseline(times, object$baseline, base_p)
  
  if (type == "survival") {
    fr_eval <- eval_laplace(b_eval$H0, object$frailty, frail_p)
    return(fr_eval$L0)
  } else if (type == "cumhaz") {
    return(b_eval$H0)
  } else {
    return(b_eval$h0)
  }
}

#' @export
predict.frailty_mle <- function(object, newdata = NULL, times = NULL,
                               type = c("survival", "cumhaz", "hazard"), ...) {
  type <- match.arg(type)
  if (is.null(times)) times <- seq(1, 100, length.out = 50)
  n_base <- get_n_base_params(object$baseline)
  n_frail <- get_n_frail_params(object$frailty)
  n_cov <- object$data$p_cov
  
  est_p <- object$coefficients
  base_p <- est_p[1:n_base]
  frail_p <- est_p[(n_base + n_cov + 1):(n_base + n_cov + n_frail)]
  
  b_eval <- eval_baseline(times, object$baseline, base_p)
  if (type == "survival") {
    fr_eval <- eval_laplace(b_eval$H0, object$frailty, frail_p)
    return(fr_eval$L0)
  } else if (type == "cumhaz") {
    return(b_eval$H0)
  } else {
    return(b_eval$h0)
  }
}
