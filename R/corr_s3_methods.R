#' S3 Methods for corrFrailty Objects
#'
#' @description
#' Printing, summary, plotting, prediction, residual extraction, and diagnostic
#' methods for fitted \code{"corrFrailty"} models.
#'
#' @name s3_methods
#' @keywords internal
NULL

#' @export
print.corrFrailty <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("
=======================================================
")
  cat("  Correlated Frailty Model for Bivariate Survival Data
")
  cat("=======================================================
")
  cat("Call:
", paste(deparse(x$call), sep = "
", collapse = "
"), "

", sep = "")
  cat("Estimation Method :", toupper(x$method), "
")
  cat("Frailty Family    :", x$frailty, "
")
  cat("Baseline 1 (T1)   :", x$baseline1, "
")
  cat("Baseline 2 (T2)   :", x$baseline2, "
")
  cat("Observations      :", x$n.obs, "pairs (Events: T1 =", x$n.events1, ", T2 =", x$n.events2, ")

")
  
  cat("Parameter Estimates:
")
  print(round(x$coefficients, digits = digits))
  
  cat("
Model Selection Criteria:
")
  crit_vec <- c(
    "Log-Likelihood" = x$logLik,
    "Deviance" = x$deviance,
    "AIC" = x$AIC,
    "BIC" = x$BIC
  )
  if (!is.na(x$DIC)) crit_vec <- c(crit_vec, "pD" = x$pD, "DIC" = x$DIC, "LPML" = x$LPML)
  print(round(crit_vec, digits = digits))
  cat("
")
  invisible(x)
}

#' @export
summary.corrFrailty <- function(object, ...) {
  res <- list(
    call = object$call,
    method = object$method,
    frailty = object$frailty,
    baseline1 = object$baseline1,
    baseline2 = object$baseline2,
    n.obs = object$n.obs,
    n.events1 = object$n.events1,
    n.events2 = object$n.events2,
    summary.table = object$summary.table,
    logLik = object$logLik,
    deviance = object$deviance,
    AIC = object$AIC,
    BIC = object$BIC,
    DIC = object$DIC,
    pD = object$pD,
    LPML = object$LPML,
    autocorrelation = object$autocorrelation
  )
  class(res) <- "summary.corrFrailty"
  res
}

#' @export
print.summary.corrFrailty <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("
=======================================================
")
  cat("  Correlated Frailty Model for Bivariate Survival Data
")
  cat("=======================================================
")
  cat("Call:
", paste(deparse(x$call), sep = "
", collapse = "
"), "

", sep = "")
  cat("Estimation Method :", toupper(x$method), "
")
  cat("Frailty Family    :", x$frailty, "
")
  cat("Baseline 1 (T1)   :", x$baseline1, "
")
  cat("Baseline 2 (T2)   :", x$baseline2, "
")
  cat("Number of Pairs   :", x$n.obs, "
")
  cat("Number of Events  : T1 =", x$n.events1, "; T2 =", x$n.events2, "

")
  
  cat("Coefficients & Posterior Summaries:
")
  print(round(x$summary.table, digits = digits))
  
  cat("
Model Selection & Comparison:
")
  cat("  Log-Likelihood :", round(x$logLik, digits), "
")
  cat("  Deviance       :", round(x$deviance, digits), "
")
  cat("  AIC            :", round(x$AIC, digits), "
")
  cat("  BIC            :", round(x$BIC, digits), "
")
  if (!is.na(x$DIC)) {
    cat("  pD (Eff. Params):", round(x$pD, digits), "
")
    cat("  DIC            :", round(x$DIC, digits), "
")
    cat("  LPML           :", round(x$LPML, digits), "
")
  }
  
  if (!is.null(x$autocorrelation)) {
    cat("
MCMC Sample Autocorrelations:
")
    print(round(x$autocorrelation, digits = digits))
  }
  cat("
")
  invisible(x)
}

#' @export
plot.corrFrailty <- function(x,
                             which = c("trace", "density", "acf", "survival", "hazard", "diagnostics"),
                             ask = prod(graphics::par("mfcol")) < length(which) && grDevices::dev.interactive(),
                             ...) {
  which <- match.arg(which, several.ok = TRUE)
  
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par))
  
  if (ask) {
    oask <- graphics::par(ask = TRUE)
    on.exit(graphics::par(oask), add = TRUE)
  }
  
  param_names <- x$param.names
  n_p <- length(param_names)
  
  # 1. Trace plots (for MCMC)
  if ("trace" %in% which) {
    if (x$method == "mcmc" && !is.null(x$all.samples)) {
      nr <- ceiling(sqrt(n_p))
      nc <- ceiling(n_p / nr)
      graphics::par(mfrow = c(nr, nc), mar = c(3, 3, 2, 1), mgp = c(1.8, 0.5, 0))
      
      for (j in seq_len(n_p)) {
        graphics::plot(x$all.samples[, j], type = "l", col = "steelblue",
                       main = paste("Trace:", param_names[j]),
                       xlab = "Iteration", ylab = "Value")
        graphics::abline(h = x$coefficients[j], col = "red", lty = 2, lwd = 1.5)
      }
    }
  }
  
  # 2. Density plots (for MCMC)
  if ("density" %in% which) {
    if (x$method == "mcmc" && !is.null(x$all.samples)) {
      nr <- ceiling(sqrt(n_p))
      nc <- ceiling(n_p / nr)
      graphics::par(mfrow = c(nr, nc), mar = c(3, 3, 2, 1), mgp = c(1.8, 0.5, 0))
      
      for (j in seq_len(n_p)) {
        d <- stats::density(x$all.samples[, j])
        graphics::plot(d, main = paste("Density:", param_names[j]),
                       col = "darkblue", lwd = 2, xlab = "Value", ylab = "Density")
        graphics::polygon(d, col = grDevices::adjustcolor("steelblue", alpha.f = 0.3), border = NA)
        graphics::abline(v = x$coefficients[j], col = "red", lty = 2, lwd = 1.5)
      }
    }
  }
  
  # 3. ACF plots (for MCMC)
  if ("acf" %in% which) {
    if (x$method == "mcmc" && !is.null(x$all.samples)) {
      nr <- ceiling(sqrt(n_p))
      nc <- ceiling(n_p / nr)
      graphics::par(mfrow = c(nr, nc), mar = c(3, 3, 2, 1), mgp = c(1.8, 0.5, 0))
      
      for (j in seq_len(n_p)) {
        acf_res <- stats::acf(x$all.samples[, j], plot = FALSE, lag.max = 30)
        graphics::plot(acf_res, main = paste("ACF:", param_names[j]))
      }
    }
  }
  
  # 4. Survival Curves
  if ("survival" %in% which) {
    graphics::par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
    max_t1 <- max(x$data$t1)
    max_t2 <- max(x$data$t2)
    grid_t1 <- seq(0.1, max_t1 * 1.2, length.out = 100)
    grid_t2 <- seq(0.1, max_t2 * 1.2, length.out = 100)
    
    unp <- x$unpack(x$coefficients)
    
    # Marginal survival curves
    H1_grid <- baseline_H_h(grid_t1, x$baseline1, unp$base1)$H
    H2_grid <- baseline_H_h(grid_t2, x$baseline2, unp$base2)$H
    
    pder1 <- log_phi_derivatives(H1_grid, rep(0, length(H1_grid)), x$frailty, unp$frailty)
    pder2 <- log_phi_derivatives(rep(0, length(H2_grid)), H2_grid, x$frailty, unp$frailty)
    
    S1_grid <- exp(pder1$phi)
    S2_grid <- exp(pder2$phi)
    
    graphics::plot(grid_t1, S1_grid, type = "l", lwd = 2, col = "blue",
                   xlab = "Time (t)", ylab = "Survival Probability S(t)",
                   main = "Estimated Marginal Survival", ylim = c(0, 1))
    graphics::lines(grid_t2, S2_grid, lwd = 2, col = "darkgreen", lty = 2)
    graphics::legend("topright", legend = c("T1 Marginal", "T2 Marginal"),
                     col = c("blue", "darkgreen"), lty = c(1, 2), lwd = 2, bty = "n")
    
    # 2D Bivariate Survival Contour Plot
    grid_len <- 30
    g1 <- seq(0.1, max_t1, length.out = grid_len)
    g2 <- seq(0.1, max_t2, length.out = grid_len)
    mesh <- expand.grid(t1 = g1, t2 = g2)
    
    h1_mesh <- baseline_H_h(mesh$t1, x$baseline1, unp$base1)$H
    h2_mesh <- baseline_H_h(mesh$t2, x$baseline2, unp$base2)$H
    
    pder_mesh <- log_phi_derivatives(h1_mesh, h2_mesh, x$frailty, unp$frailty)
    S_mesh <- matrix(exp(pder_mesh$phi), nrow = grid_len, ncol = grid_len)
    
    graphics::contour(g1, g2, S_mesh, nlevels = 10, col = "purple",
                      xlab = "Time t1", ylab = "Time t2",
                      main = "Bivariate Survival Contours S(t1, t2)")
  }
  
  # 5. Baseline Hazard Curves
  if ("hazard" %in% which) {
    graphics::par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
    max_t1 <- max(x$data$t1)
    max_t2 <- max(x$data$t2)
    grid_t1 <- seq(0.1, max_t1 * 1.2, length.out = 100)
    grid_t2 <- seq(0.1, max_t2 * 1.2, length.out = 100)
    
    unp <- x$unpack(x$coefficients)
    h1_grid <- baseline_H_h(grid_t1, x$baseline1, unp$base1)$h
    h2_grid <- baseline_H_h(grid_t2, x$baseline2, unp$base2)$h
    
    graphics::plot(grid_t1, h1_grid, type = "l", lwd = 2, col = "firebrick",
                   xlab = "Time (t1)", ylab = "Hazard Rate h_01(t1)",
                   main = paste("Baseline Hazard T1 (", x$baseline1, ")", sep = ""))
    graphics::plot(grid_t2, h2_grid, type = "l", lwd = 2, col = "darkorange",
                   xlab = "Time (t2)", ylab = "Hazard Rate h_02(t2)",
                   main = paste("Baseline Hazard T2 (", x$baseline2, ")", sep = ""))
  }
  
  # 6. Diagnostics (Geweke Z plot)
  if ("diagnostics" %in% which && x$method == "mcmc") {
    graphics::par(mfrow = c(1, 1), mar = c(5, 5, 3, 2))
    gz <- x$summary.table$Geweke.Z
    graphics::plot(seq_len(n_p), gz, type = "p", pch = 19, col = "darkblue", cex = 1.3,
                   xlab = "Parameter", ylab = "Geweke Z-score",
                   main = "Geweke Convergence Diagnostic Z-Scores",
                   xaxt = "n", ylim = c(min(-2.5, min(gz) - 0.5), max(2.5, max(gz) + 0.5)))
    graphics::axis(1, at = seq_len(n_p), labels = param_names, las = 2)
    graphics::abline(h = 0, col = "gray50", lty = 1)
    graphics::abline(h = c(-1.96, 1.96), col = "red", lty = 2, lwd = 1.5)
  }
  
  invisible(x)
}

#' @export
predict.corrFrailty <- function(object, newdata = NULL, times = NULL,
                                type = c("survival", "hazard", "cumulative_hazard"), ...) {
  type <- match.arg(type)
  unp <- object$unpack(object$coefficients)
  
  if (is.null(times)) {
    times <- seq(0.1, max(object$data$t1, object$data$t2), length.out = 50)
  }
  
  # Linear predictor for newdata
  eta <- if (!is.null(newdata) && length(unp$beta) > 0) {
    X_new <- stats::model.matrix(object$formula[-2], data = newdata)
    if ("(Intercept)" %in% colnames(X_new)) X_new <- X_new[, colnames(X_new) != "(Intercept)", drop = FALSE]
    as.vector(exp(X_new %*% unp$beta))
  } else {
    rep(1, length(times))
  }
  
  res1 <- baseline_H_h(times, object$baseline1, unp$base1)
  res2 <- baseline_H_h(times, object$baseline2, unp$base2)
  
  if (type == "hazard") {
    return(data.frame(time = times, hazard1 = res1$h * mean(eta), hazard2 = res2$h * mean(eta)))
  } else if (type == "cumulative_hazard") {
    return(data.frame(time = times, cumhaz1 = res1$H * mean(eta), cumhaz2 = res2$H * mean(eta)))
  } else {
    # Marginal survival
    s1 <- res1$H * mean(eta)
    s2 <- res2$H * mean(eta)
    pder1 <- log_phi_derivatives(s1, rep(0, length(s1)), object$frailty, unp$frailty)
    pder2 <- log_phi_derivatives(rep(0, length(s2)), s2, object$frailty, unp$frailty)
    return(data.frame(time = times, survival1 = exp(pder1$phi), survival2 = exp(pder2$phi)))
  }
}

#' @export
residuals.corrFrailty <- function(object, type = c("cox-snell", "martingale"), ...) {
  type <- match.arg(type)
  unp <- object$unpack(object$coefficients)
  
  t1 <- object$data$t1; d1 <- object$data$d1
  t2 <- object$data$t2; d2 <- object$data$d2
  X <- object$data$X
  
  eta <- if (!is.null(X) && length(unp$beta) > 0) {
    as.vector(exp(X %*% unp$beta))
  } else {
    rep(1, length(t1))
  }
  
  H1 <- baseline_H_h(t1, object$baseline1, unp$base1)$H * eta
  H2 <- baseline_H_h(t2, object$baseline2, unp$base2)$H * eta
  
  if (type == "cox-snell") {
    return(data.frame(cs1 = H1, cs2 = H2))
  } else {
    return(data.frame(mart1 = d1 - H1, mart2 = d2 - H2))
  }
}

#' @export
coef.corrFrailty <- function(object, ...) {
  object$coefficients
}

#' @export
vcov.corrFrailty <- function(object, ...) {
  if (object$method == "mcmc") {
    stats::var(object$all.samples)
  } else {
    object$vcov
  }
}

#' @export
logLik.corrFrailty <- function(object, ...) {
  val <- object$logLik
  attr(val, "df") <- object$n.params
  attr(val, "nobs") <- object$n.obs
  class(val) <- "logLik"
  val
}

#' @export
AIC.corrFrailty <- function(object, ...) {
  object$AIC
}

#' @export
BIC.corrFrailty <- function(object, ...) {
  object$BIC
}
