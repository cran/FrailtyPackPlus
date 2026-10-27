#' Diagnostic and Visual Plots for One-Shot Frailty Models
#'
#' @description
#' Generates publication-ready diagnostic plots for fitted \code{"oneshot_kofm"} models,
#' including \eqn{k}-out-of-\eqn{M} system reliability curves, component marginal survival
#' curves, observed vs. fitted cell proportions, and standardized Pearson residuals.
#'
#' @param x A fitted model object of class \code{"oneshot_kofm"}.
#' @param which A character string or integer specifying the plot type:
#'   \code{"all"} (default), \code{"reliability"} (or 1), \code{"marginal"} (or 2),
#'   \code{"fit"} (or 3), or \code{"residuals"} (or 4).
#' @param stress A numeric scalar specifying the operating stress condition for reliability curves.
#'   If \code{NULL} (default), uses the minimum stress level in the dataset.
#' @param t_max A positive numeric scalar specifying the maximum mission time for curves.
#'   If \code{NULL} (default), uses \eqn{1.5 \times \max(\tau_i)}.
#' @param col Optional vector of colors for curves.
#' @param lty Optional vector of line types.
#' @param lwd Line width (default 2).
#' @param ... Additional graphical parameters passed to plotting functions.
#'
#' @return No return value, called for plotting side effects.
#'
#' @export
plot.oneshot_kofm <- function(x, which = c("all", "reliability", "marginal", "fit", "residuals"),
                              stress = NULL, t_max = NULL, col = NULL, lty = NULL, lwd = 2, ...) {
  if (!inherits(x, "oneshot_kofm")) {
    stop("'x' must be a fitted object of class 'oneshot_kofm'.")
  }
  
  if (is.numeric(which)) {
    type_map <- c("reliability", "marginal", "fit", "residuals")
    if (which >= 1 && which <= 4) which <- type_map[which] else which <- "all"
  } else {
    which <- match.arg(which)
  }
  
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par))
  
  tau_max <- max(x$data$time, na.rm = TRUE)
  if (is.null(t_max) || t_max <= 0) t_max <- 1.5 * tau_max
  
  if (is.null(stress)) {
    stress <- if ("stress" %in% colnames(x$data)) min(x$data$stress, na.rm = TRUE) else 0
  }
  
  t_seq <- seq(0, t_max, length.out = 150)
  M <- x$M
  
  if (which == "all") {
    graphics::par(mfrow = c(2, 2), mar = c(4.2, 4.2, 2.5, 1.5))
  }
  
  # 1. System Reliability Curves
  if (which %in% c("all", "reliability")) {
    rel_df <- reliability_k_out_of_m(x, t = t_seq, stress = stress)
    
    if (is.null(col)) {
      pal <- grDevices::palette.colors(n = M, palette = "Okabe-Ito")
    } else {
      pal <- rep(col, length.out = M)
    }
    line_types <- if (is.null(lty)) seq_len(M) else rep(lty, length.out = M)
    
    graphics::plot(0, 0, type = "n", xlim = c(0, t_max), ylim = c(0, 1),
                   xlab = "Time (t)", ylab = expression(R[k](t)),
                   main = paste0("k-out-of-", M, " System Reliability (s = ", round(stress, 1), ")"), ...)
    graphics::grid(col = "gray90")
    
    for (k_idx in seq_len(M)) {
      sub_df <- rel_df[rel_df$k == k_idx, ]
      graphics::lines(sub_df$time, sub_df$reliability, col = pal[k_idx],
                      lty = line_types[k_idx], lwd = lwd)
    }
    
    legend_labels <- vapply(seq_len(M), function(k_i) {
      if (k_i == 1) "k=1 (Parallel)"
      else if (k_i == M) paste0("k=", M, " (Series)")
      else paste0("k=", k_i)
    }, character(1))
    
    graphics::legend("topright", legend = legend_labels, col = pal,
                     lty = line_types, lwd = lwd, bty = "n", cex = 0.85)
  }
  
  # 2. Marginal Component Survival Curves
  if (which %in% c("all", "marginal")) {
    if (is.null(col)) {
      pal <- grDevices::palette.colors(n = M, palette = "Okabe-Ito")
    } else {
      pal <- rep(col, length.out = M)
    }
    line_types <- if (is.null(lty)) seq_len(M) else rep(lty, length.out = M)
    
    moms <- component_moments(x, stress = stress)
    beta_val <- if (x$frailty) x$coefficients["beta"] else 0
    
    graphics::plot(0, 0, type = "n", xlim = c(0, t_max), ylim = c(0, 1),
                   xlab = "Time (t)", ylab = "Marginal Survival S(t)",
                   main = paste0("Component Survival Curves (s = ", round(stress, 1), ")"), ...)
    graphics::grid(col = "gray90")
    
    for (m in seq_len(M)) {
      lam_m <- moms$lambda[m]
      surv_m <- if (beta_val > 1e-10) {
        (1 + beta_val * lam_m * t_seq)^(-1 / beta_val)
      } else {
        exp(-lam_m * t_seq)
      }
      graphics::lines(t_seq, surv_m, col = pal[m], lty = line_types[m], lwd = lwd)
    }
    
    graphics::legend("topright", legend = x$components, col = pal,
                     lty = line_types, lwd = lwd, bty = "n", cex = 0.85)
  }
  
  # 3. Model Calibration: Observed vs Fitted Proportions
  if (which %in% c("all", "fit")) {
    obs_counts <- x$data[, colnames(x$data) %in% colnames(x$fitted_counts), drop = FALSE]
    if (ncol(obs_counts) == 0) {
      obs_counts <- x$fitted_counts + x$residuals * sqrt(pmax(x$fitted_counts, 1e-8))
    }
    K_i <- rowSums(obs_counts)
    obs_probs <- as.matrix(obs_counts) / ifelse(K_i > 0, K_i, 1)
    fit_probs <- x$fitted_probs
    
    graphics::plot(as.vector(fit_probs), as.vector(obs_probs),
                   pch = 19, col = "royalblue",
                   xlab = "Fitted Subset Probabilities",
                   ylab = "Observed Subset Proportions",
                   main = "Model Fit: Observed vs. Fitted", ...)
    graphics::abline(0, 1, col = "firebrick", lty = 2, lwd = 1.5)
    graphics::grid(col = "gray90")
  }
  
  # 4. Standardized Pearson Residuals
  if (which %in% c("all", "residuals")) {
    res_vec <- as.vector(x$residuals)
    graphics::plot(res_vec, pch = 19, col = "darkorange",
                   xlab = "Cell Index",
                   ylab = "Standardized Pearson Residual",
                   main = "Residual Diagnostics", ...)
    graphics::abline(h = 0, col = "gray40", lty = 1)
    graphics::abline(h = c(-2, 2), col = "firebrick", lty = 2)
    graphics::grid(col = "gray90")
  }
  
  invisible(x)
}


#' Diagnostic and Reliability Plots for One-Shot Device Testing
#'
#' @description
#' Provides diagnostic visualization and reliability plotting for fitted
#' \code{"oneshot_frailty"} model objects.
#'
#' @param x A fitted object of class \code{"oneshot_frailty"}.
#' @param which A character string or integer specifying the plot type:
#'   \code{"all"} (default), \code{"reliability"} (1), \code{"survival"} (2),
#'   \code{"residuals"} (3), or \code{"fitted"} (4).
#' @param stress Optional numeric vector or value of stress covariate for reliability evaluation.
#' @param t_max Optional maximum time for reliability curves (defaults to \eqn{1.5 x \max(\tau)}).
#' @param ask Logical: \code{TRUE} to prompt before new plot pages in interactive mode.
#' @param ... Additional graphical parameters passed to \code{plot}.
#'
#' @return Invisible \code{NULL}. Produces base R graphical plots.
#'
#' @details
#' The available plot types are:
#' \describe{
#'   \item{\code{"reliability"}}{System reliability curves \eqn{R_k(t)} for \eqn{k=1} (parallel),
#'     \eqn{k=M} (series), and intermediate \eqn{k}-out-of-\eqn{M} structures.}
#'   \item{\code{"survival"}}{Marginal component survival curves \eqn{S_m(t)}.}
#'   \item{\code{"residuals"}}{Standardized Pearson residuals vs. test group index.}
#'   \item{\code{"fitted"}}{Observed vs. model-fitted device counts for each malfunction subset.}
#' }
#'
#' @references
#' Yu, C., So, H. Y., and Balakrishnan, N. (2025). \emph{Inference for one-shot devices
#' with Weibull component lifetimes under gamma frailty model}. In: Stochastic Modeling
#' and Statistical Methods, Elsevier Inc., Chapter 1, pp. 1--34.
#' \doi{10.1016/B978-0-44-331694-4.00006-2}.
#'
#' @examples
#' data(classH)
#' start_par <- c(-3.997, -0.009, -14.412, 0.003, 0.004, -0.001, -5.493e-4, -0.005, 0.082)
#' fit_w <- oneshot_frailty(
#'   data = classH,
#'   time = "time",
#'   stress = "stress",
#'   counts = c("n_none", "n_T", "n_G", "n_both"),
#'   baseline = "weibull",
#'   frailty = TRUE,
#'   start = start_par
#' )
#' plot(fit_w, which = "reliability", stress = 400)
#' @export
plot.oneshot_frailty <- function(x,
                                 which = c("all", "reliability", "survival", "residuals", "fitted"),
                                 stress = NULL,
                                 t_max = NULL,
                                 ask = prod(graphics::par("mfcol")) < 2 && grDevices::dev.interactive(),
                                 ...) {
  if (is.numeric(which)) {
    which_map <- c("reliability", "survival", "residuals", "fitted")
    if (which >= 1 && which <= 4) which <- which_map[which]
    else which <- "all"
  } else {
    which <- match.arg(which)
  }
  
  if (ask) {
    oask <- grDevices::devAskNewPage(TRUE)
    on.exit(grDevices::devAskNewPage(oask))
  }
  
  if (is.null(t_max)) {
    t_max <- max(x$tau) * 1.5
  }
  t_seq <- seq(0, t_max, length.out = 100)
  
  # Plot 1: Reliability curves
  if (which %in% c("all", "reliability")) {
    rel_df <- reliability_k_out_of_m(x, t = t_seq, stress = stress)
    k_vals <- unique(rel_df$k)
    colors <- grDevices::palette.colors(n = max(length(k_vals), 2), palette = "Dark 2")
    
    graphics::plot(t_seq, rel_df$reliability[rel_df$k == k_vals[1]], type = "n",
                   ylim = c(0, 1), xlab = "Mission Time (t)", ylab = "System Reliability R_k(t)",
                   main = sprintf("k-out-of-%d System Reliability (%s Baseline)", x$M,
                                  if (x$baseline == "weibull") "Weibull" else "Exponential"), ...)
    graphics::grid()
    
    for (i in seq_along(k_vals)) {
      kv <- k_vals[i]
      sub_d <- rel_df[rel_df$k == kv, ]
      graphics::lines(sub_d$time, sub_d$reliability, col = colors[i], lwd = 2, lty = i)
    }
    
    leg_labels <- vapply(k_vals, function(k) {
      if (k == 1) "Parallel (k=1)"
      else if (k == x$M) sprintf("Series (k=%d)", x$M)
      else sprintf("k=%d", k)
    }, character(1))
    
    graphics::legend("topright", legend = leg_labels, col = colors[seq_along(k_vals)],
                     lwd = 2, lty = seq_along(k_vals), bty = "n")
  }
  
  # Plot 2: Marginal survival curves
  if (which %in% c("all", "survival")) {
    p_cov <- ncol(x$X_cov)
    if (is.null(stress)) {
      X_eval <- matrix(colMeans(x$X_cov), nrow = 1)
    } else {
      if (length(stress) == p_cov - 1) X_eval <- matrix(c(1, stress), nrow = 1)
      else X_eval <- matrix(stress, nrow = 1)
    }
    
    colors <- grDevices::palette.colors(n = max(x$M, 2), palette = "Dark 2")
    graphics::plot(t_seq, rep(1, length(t_seq)), type = "n", ylim = c(0, 1),
                   xlab = "Time (t)", ylab = "Marginal Survival Probability S_m(t)",
                   main = "Marginal Component Survival Functions", ...)
    graphics::grid()
    
    for (m in seq_len(x$M)) {
      s_vec <- vapply(t_seq, function(tt) {
        pr <- compute_probabilities(x$coefficients, X_eval, tt, x$M, x$baseline, x$frailty)
        lam <- pr$lambda_mat[1, m]
        alp <- pr$alpha_mat[1, m]
        beta_v <- pr$beta
        h0 <- (lam * tt)^alp
        if (beta_v > 1e-10) (1 + beta_v * h0)^(-1 / beta_v)
        else exp(-h0)
      }, numeric(1))
      graphics::lines(t_seq, s_vec, col = colors[m], lwd = 2, lty = m)
    }
    graphics::legend("topright", legend = x$components, col = colors[seq_len(x$M)],
                     lwd = 2, lty = seq_len(x$M), bty = "n")
  }
  
  # Plot 3: Pearson Residuals
  if (which %in% c("all", "residuals")) {
    res_mat <- x$residuals
    n_groups <- nrow(res_mat)
    n_subsets <- ncol(res_mat)
    
    graphics::plot(seq_len(n_groups), rep(0, n_groups), type = "n",
                   ylim = c(min(-3, min(res_mat)), max(3, max(res_mat))),
                   xlab = "Test Group Index", ylab = "Standardized Pearson Residual",
                   main = "Pearson Residuals by Inspection Group", ...)
    graphics::grid()
    graphics::abline(h = c(-2, 0, 2), col = c("red", "black", "red"), lty = c(2, 1, 2))
    
    colors <- grDevices::palette.colors(n = max(n_subsets, 2), palette = "Dark 2")
    for (j in seq_len(n_subsets)) {
      graphics::points(seq_len(n_groups), res_mat[, j], col = colors[j], pch = 15 + j, cex = 1.2)
    }
    graphics::legend("topright", legend = colnames(res_mat), col = colors[seq_len(n_subsets)],
                     pch = 15 + seq_len(n_subsets), bty = "n")
  }
  
  # Plot 4: Observed vs. Fitted Counts
  if (which %in% c("all", "fitted")) {
    obs_counts <- as.vector(x$count_mat)
    fit_counts <- as.vector(x$fitted_counts)
    max_c <- max(c(obs_counts, fit_counts), na.rm = TRUE)
    
    graphics::plot(fit_counts, obs_counts, pch = 19, col = "navy",
                   xlim = c(0, max_c * 1.1), ylim = c(0, max_c * 1.1),
                   xlab = "Fitted Counts", ylab = "Observed Counts",
                   main = "Observed vs. Expected Device Counts", ...)
    graphics::grid()
    graphics::abline(a = 0, b = 1, col = "red", lty = 2, lwd = 2)
  }
  
  invisible(NULL)
}
