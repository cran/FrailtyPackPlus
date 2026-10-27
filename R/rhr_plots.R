# S3 Plot Methods for Shared RHR Frailty Models

#' Plot Diagnostics for Bayesian Shared RHR Frailty Model
#'
#' Produces graphical diagnostics for fitted Bayesian shared frailty models,
#' including trace plots, posterior density plots with HPD intervals, autocorrelation decay,
#' Kolmogorov-Smirnov goodness-of-fit plots, conditional predictive ordinate (CPO) plots,
#' and reversed hazard rate (RHR) curves.
#'
#' @param x An object of class \code{"rhr_bayes_fit"}.
#' @param type Character vector or string specifying plot type. One or more of:
#'   \code{"trace"}, \code{"density"}, \code{"autocorr"}, \code{"ks"},
#'   \code{"cpo"}, \code{"rhr"}. Default is \code{c("trace", "density")}.
#' @param params Optional character vector of parameter names to plot.
#' @param ask Logical; if \code{TRUE}, user is prompted before each new page of plots. Default is \code{FALSE}.
#' @param ... Additional arguments passed to base plotting functions.
#'
#' @return Invisible \code{NULL}.
#' @importFrom graphics plot par abline lines legend points
#' @export
#'
#' @examples
#' data(austrian_twin)
#' fit <- fit_rhr_frailty_bayes(cbind(t1, t2, d1, d2) ~ age, data = austrian_twin[1:30, ],
#'                               frailty = "WL", baseline = "GE", n_iter = 40, burn_in = 10)
#' plot(fit, type = "trace", params = "theta")
plot.rhr_bayes_fit <- function(x, type = c("trace", "density", "autocorr", "ks", "cpo", "rhr"),
                               params = NULL, ask = FALSE, ...) {
  if (missing(type)) type <- c("trace", "density")
  draws <- x$draws
  all_names <- colnames(draws)
  if (is.null(params)) params <- all_names
  params <- intersect(params, all_names)
  if (length(params) == 0) params <- all_names
  
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par))
  
  if (ask) graphics::par(ask = TRUE)
  
  # 1. Trace Plots
  if ("trace" %in% type) {
    n_p <- length(params)
    n_col <- min(3, n_p)
    n_row <- ceiling(n_p / n_col)
    graphics::par(mfrow = c(n_row, n_col), mar = c(3.5, 3.5, 2.5, 1))
    
    for (p in params) {
      y <- draws[, p]
      graphics::plot(y, type = "l", col = "steelblue", xlab = "Iteration", ylab = p,
                     main = paste0("Trace of ", p), ...)
      graphics::abline(h = mean(y), col = "firebrick", lty = 2, lwd = 1.5)
    }
  }
  
  # 2. Posterior Density Plots
  if ("density" %in% type) {
    n_p <- length(params)
    n_col <- min(3, n_p)
    n_row <- ceiling(n_p / n_col)
    graphics::par(mfrow = c(n_row, n_col), mar = c(3.5, 3.5, 2.5, 1))
    
    for (p in params) {
      y <- draws[, p]
      d_est <- stats::density(y)
      hpd <- hpd_interval(y, 0.95)
      
      graphics::plot(d_est, col = "navy", lwd = 2, xlab = p, ylab = "Posterior Density",
                     main = paste0("Density of ", p), ...)
      graphics::abline(v = mean(y), col = "firebrick", lty = 2, lwd = 1.5)
      graphics::abline(v = hpd, col = "darkgreen", lty = 3, lwd = 1.5)
      graphics::legend("topright", legend = c("Mean", "95% HPD"), col = c("firebrick", "darkgreen"),
                       lty = c(2, 3), bty = "n", cex = 0.8)
    }
  }
  
  # 3. Autocorrelation Plots
  if ("autocorr" %in% type) {
    n_p <- length(params)
    n_col <- min(3, n_p)
    n_row <- ceiling(n_p / n_col)
    graphics::par(mfrow = c(n_row, n_col), mar = c(3.5, 3.5, 2.5, 1))
    
    for (p in params) {
      stats::acf(draws[, p], main = paste0("ACF of ", p), col = "darkslategray", ...)
    }
  }
  
  # 4. Kolmogorov-Smirnov / CDF Goodness-of-Fit
  if ("ks" %in% type) {
    graphics::par(mfrow = c(1, 2), mar = c(4, 4, 2.5, 1))
    t1_obs <- x$data$t1[x$data$d1 == 1]
    t2_obs <- x$data$t2[x$data$d2 == 1]
    
    base_dist <- x$model_info$baseline
    p_self <- x$estimates$SELF
    
    b_names <- switch(base_dist,
      "GE" = c("zeta", "delta"),
      "EG" = c("zeta", "delta"),
      "GIR" = c("zeta", "delta"),
      "MIW" = c("zeta", "delta", "xi"),
      "GLL-I" = c("lambda", "alpha", "gamma"),
      "GLL-II" = c("lambda", "alpha", "gamma"),
      "GIE" = c("zeta", "delta"),
      "Weibull" = c("zeta", "delta"),
      "Exponential" = c("delta"),
      "Log-Logistic" = c("lambda", "gamma")
    )
    
    p1_est <- p_self[paste0(b_names, "1")]
    names(p1_est) <- b_names
    p2_est <- p_self[paste0(b_names, "2")]
    names(p2_est) <- b_names
    
    if (length(t1_obs) > 0) {
      grid1 <- seq(min(t1_obs), max(t1_obs), length.out = 100)
      emp_cdf1 <- stats::ecdf(t1_obs)
      fit_cdf1 <- p_rhr_baseline(grid1, base_dist, p1_est)
      
      graphics::plot(emp_cdf1, col = "gray40", lwd = 1.5, xlab = "T1", ylab = "Cumulative Probability",
                     main = "K-S / CDF Fit for T1")
      graphics::lines(grid1, fit_cdf1, col = "firebrick", lwd = 2)
      graphics::legend("bottomright", legend = c("Empirical", "Model Fit"),
                       col = c("gray40", "firebrick"), lty = c(1, 1), lwd = c(1.5, 2), bty = "n")
    }
    
    if (length(t2_obs) > 0) {
      grid2 <- seq(min(t2_obs), max(t2_obs), length.out = 100)
      emp_cdf2 <- stats::ecdf(t2_obs)
      fit_cdf2 <- p_rhr_baseline(grid2, base_dist, p2_est)
      
      graphics::plot(emp_cdf2, col = "gray40", lwd = 1.5, xlab = "T2", ylab = "Cumulative Probability",
                     main = "K-S / CDF Fit for T2")
      graphics::lines(grid2, fit_cdf2, col = "navy", lwd = 2)
      graphics::legend("bottomright", legend = c("Empirical", "Model Fit"),
                       col = c("gray40", "navy"), lty = c(1, 1), lwd = c(1.5, 2), bty = "n")
    }
  }
  
  # 5. CPO Diagnostic Plot
  if ("cpo" %in% type && !is.null(x$cpo)) {
    graphics::par(mfrow = c(1, 1), mar = c(4, 4, 3, 1))
    cpo_vals <- x$cpo
    graphics::plot(cpo_vals, type = "h", col = "steelblue", xlab = "Observation Index (j)",
                   ylab = "Conditional Predictive Ordinate (CPO)",
                   main = "CPO Index Plot for Model Diagnostic", ...)
    graphics::points(cpo_vals, pch = 19, cex = 0.6, col = "navy")
    graphics::abline(h = mean(cpo_vals), col = "firebrick", lty = 2)
  }
  
  # 6. Reversed Hazard Rate (RHR) Curves
  if ("rhr" %in% type) {
    graphics::par(mfrow = c(1, 2), mar = c(4, 4, 2.5, 1))
    t1_obs <- x$data$t1
    t2_obs <- x$data$t2
    base_dist <- x$model_info$baseline
    p_self <- x$estimates$SELF
    
    b_names <- switch(base_dist,
      "GE" = c("zeta", "delta"),
      "EG" = c("zeta", "delta"),
      "GIR" = c("zeta", "delta"),
      "MIW" = c("zeta", "delta", "xi"),
      "GLL-I" = c("lambda", "alpha", "gamma"),
      "GLL-II" = c("lambda", "alpha", "gamma"),
      "GIE" = c("zeta", "delta"),
      "Weibull" = c("zeta", "delta"),
      "Exponential" = c("delta"),
      "Log-Logistic" = c("lambda", "gamma")
    )
    
    p1_est <- p_self[paste0(b_names, "1")]
    names(p1_est) <- b_names
    p2_est <- p_self[paste0(b_names, "2")]
    names(p2_est) <- b_names
    
    g1 <- seq(min(t1_obs), max(t1_obs), length.out = 100)
    rhr1 <- rhr_baseline(g1, base_dist, p1_est)
    graphics::plot(g1, rhr1, type = "l", col = "firebrick", lwd = 2,
                   xlab = "Time (T1)", ylab = "Reversed Hazard Rate r01(t)",
                   main = "Estimated RHR for T1")
    
    g2 <- seq(min(t2_obs), max(t2_obs), length.out = 100)
    rhr2 <- rhr_baseline(g2, base_dist, p2_est)
    graphics::plot(g2, rhr2, type = "l", col = "navy", lwd = 2,
                   xlab = "Time (T2)", ylab = "Reversed Hazard Rate r02(t)",
                   main = "Estimated RHR for T2")
  }
  
  invisible(NULL)
}

#' Plot Diagnostics for Maximum Likelihood Shared RHR Frailty Model
#'
#' Produces goodness-of-fit and reversed hazard curves for an MLE shared frailty model.
#'
#' @param x An object of class \code{"rhr_mle_fit"}.
#' @param type Character string indicating plot type: \code{"ks"} or \code{"rhr"}.
#' @param ... Additional arguments passed to plotting functions.
#'
#' @return Invisible \code{NULL}.
#' @importFrom graphics plot par lines
#' @export
#'
#' @examples
#' data(austrian_twin)
#' fit <- fit_rhr_frailty_mle(cbind(t1, t2, d1, d2) ~ age, data = austrian_twin[1:30, ],
#'                            frailty = "WL", baseline = "GE")
#' plot(fit, type = "ks")
plot.rhr_mle_fit <- function(x, type = c("ks", "rhr"), ...) {
  type <- match.arg(type)
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par))
  
  t1_obs <- x$data$t1[x$data$d1 == 1]
  t2_obs <- x$data$t2[x$data$d2 == 1]
  base_dist <- x$model_info$baseline
  coefs <- x$coefficients
  
  b_names <- switch(base_dist,
    "GE" = c("zeta", "delta"),
    "EG" = c("zeta", "delta"),
    "GIR" = c("zeta", "delta"),
    "MIW" = c("zeta", "delta", "xi"),
    "GLL-I" = c("lambda", "alpha", "gamma"),
    "GLL-II" = c("lambda", "alpha", "gamma"),
    "GIE" = c("zeta", "delta"),
    "Weibull" = c("zeta", "delta"),
    "Exponential" = c("delta"),
    "Log-Logistic" = c("lambda", "gamma")
  )
  
  p1_est <- coefs[paste0(b_names, "1")]
  names(p1_est) <- b_names
  p2_est <- coefs[paste0(b_names, "2")]
  names(p2_est) <- b_names
  
  if (type == "ks") {
    graphics::par(mfrow = c(1, 2), mar = c(4, 4, 2.5, 1))
    grid1 <- seq(min(t1_obs), max(t1_obs), length.out = 100)
    graphics::plot(stats::ecdf(t1_obs), col = "gray40", lwd = 1.5, xlab = "T1", ylab = "CDF", main = "K-S Fit for T1")
    graphics::lines(grid1, p_rhr_baseline(grid1, base_dist, p1_est), col = "firebrick", lwd = 2)
    
    grid2 <- seq(min(t2_obs), max(t2_obs), length.out = 100)
    graphics::plot(stats::ecdf(t2_obs), col = "gray40", lwd = 1.5, xlab = "T2", ylab = "CDF", main = "K-S Fit for T2")
    graphics::lines(grid2, p_rhr_baseline(grid2, base_dist, p2_est), col = "navy", lwd = 2)
  } else if (type == "rhr") {
    graphics::par(mfrow = c(1, 2), mar = c(4, 4, 2.5, 1))
    g1 <- seq(min(t1_obs), max(t1_obs), length.out = 100)
    graphics::plot(g1, rhr_baseline(g1, base_dist, p1_est), type = "l", col = "firebrick", lwd = 2,
                   xlab = "T1", ylab = "r01(t)", main = "Estimated RHR for T1")
    g2 <- seq(min(t2_obs), max(t2_obs), length.out = 100)
    graphics::plot(g2, rhr_baseline(g2, base_dist, p2_est), type = "l", col = "navy", lwd = 2,
                   xlab = "T2", ylab = "r02(t)", main = "Estimated RHR for T2")
  }
  
  invisible(NULL)
}
