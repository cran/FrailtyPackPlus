#' Plot Diagnostics and Survival Curves for ziprom_fit Objects
#'
#' Provides rich visual diagnostics including fitted population and susceptible
#' survival curves, multi-level significance forecast plots (90\%, 95\%, 99\%),
#' normal QQ-plots of randomized quantile residuals with simulated envelopes,
#' residual index plots, hazard/density plots, and activation scheme comparisons.
#'
#' @param x An object of class \code{"ziprom_fit"}.
#' @param which Character string specifying the plot type:
#'   \itemize{
#'     \item \code{"survival"} (default): Fitted population survival curve \eqn{S_p(t)} with Kaplan-Meier overlay.
#'     \item \code{"forecast"}: Forecast survival curve with 90\%, 95\%, and 99\% confidence bands.
#'     \item \code{"qq"}: Normal QQ-plot of randomized quantile residuals with 95\% confidence envelope.
#'     \item \code{"residuals"}: Residuals versus observation index or fitted values.
#'     \item \code{"hazard"}: Baseline hazard rate curve.
#'     \item \code{"density"}: Population density curve.
#'     \item \code{"activation_compare"}: Comparative survival curves for Minimum, Maximum, and Random schemes.
#'   }
#' @param newdata Optional data frame containing covariate values for plotting.
#' @param times Optional numeric vector of time points for plotting curves.
#' @param level Numeric vector of confidence levels for forecast bands. Default is \code{c(0.90, 0.95, 0.99)}.
#' @param col Vector of colors used for curves and confidence bands.
#' @param lty Vector of line types.
#' @param lwd Line width (default \code{lwd = 2}).
#' @param main Custom title for the plot.
#' @param xlab Custom x-axis label.
#' @param ylab Custom y-axis label.
#' @param ask Logical; if \code{TRUE}, the user is asked before each plot. Default is \code{FALSE}.
#' @param ... Additional arguments passed to base graphics functions.
#'
#' @return Invisible \code{NULL}. The function is called for its graphical side effects.
#'
#' @references
#' Assuncao, D., Ramos, P. L., Tomazella, V., Leao, J., Osses, L., and Rodrigues, A. (2026).
#' Zero-Inflated Promotion Cure-Rate Models With Frailty and Flexible Activation Schemes.
#' \emph{Biometrical Journal}, 68:e70162. \doi{10.1002/bimj.70162}
#'
#' @importFrom stats predict
#' @export
#' @examples
#' set.seed(42)
#' sim_dat <- sim_ziprom(n = 60, activation = "random")
#' fit <- ziprom_fit(survival::Surv(time, status) ~ x, data = sim_dat, activation = "random")
#' plot(fit, which = "survival")
#' plot(fit, which = "forecast")
#' plot(fit, which = "qq")
plot.ziprom_fit <- function(x, which = c("survival", "forecast", "qq", "residuals",
                                        "hazard", "density", "activation_compare"),
                            newdata = NULL, times = NULL,
                            level = c(0.90, 0.95, 0.99),
                            col = NULL, lty = NULL, lwd = 2,
                            main = NULL, xlab = NULL, ylab = NULL,
                            ask = FALSE, ...) {
  which <- match.arg(which)

  max_t <- if (length(x$time) > 0) max(x$time) else 10
  if (is.null(times)) {
    times <- seq(0, max_t, length.out = 150)
  }
  times <- sort(as.numeric(times))

  if (which == "survival") {
    pred_s <- predict(x, newdata = newdata, type = "survival", times = times)
    s_curve <- if (is.matrix(pred_s)) colMeans(pred_s) else pred_s

    if (is.null(col)) col <- "navy"
    if (is.null(xlab)) xlab <- "Time"
    if (is.null(ylab)) ylab <- "Population Survival Probability Sp(t)"
    if (is.null(main)) main <- paste("Fitted Survival Curve -", toupper(x$activation), "Scheme")

    graphics::plot(times, s_curve, type = "l", col = col, lwd = lwd,
                   ylim = c(0, 1), xlab = xlab, ylab = ylab, main = main, ...)
    graphics::grid(col = "gray90", lty = 1)
    graphics::lines(times, s_curve, col = col, lwd = lwd)

    # Overlay Kaplan-Meier if raw data is available
    if (length(x$time) > 0 && length(x$status) > 0) {
      km_fit <- survival::survfit(survival::Surv(x$time, x$status) ~ 1)
      graphics::lines(km_fit, col = "gray40", lty = 2, lwd = 1.5, mark.time = FALSE)
      graphics::legend("topright",
                       legend = c(paste("Model (", x$activation, ")", sep = ""), "Kaplan-Meier"),
                       col = c(col, "gray40"), lty = c(1, 2), lwd = c(lwd, 1.5),
                       bg = "white", bty = "o")
    }
  } else if (which == "forecast") {
    fc_df <- forecast_ziprom(x, newdata = newdata, times = times, level = level)

    if (is.null(xlab)) xlab <- "Forecast Horizon (Time)"
    if (is.null(ylab)) ylab <- "Forecast Survival Probability"
    if (is.null(main)) main <- "Multi-Level Significance Survival Forecast"

    # Plot base curve
    graphics::plot(fc_df$time, fc_df$survival, type = "n",
                   ylim = c(0, 1), xlab = xlab, ylab = ylab, main = main, ...)
    graphics::grid(col = "gray90", lty = 1)

    # Shaded confidence bands (99%, 95%, 90%)
    graphics::polygon(c(fc_df$time, rev(fc_df$time)),
                      c(fc_df$upper_99, rev(fc_df$lower_99)),
                      col = grDevices::adjustcolor("steelblue", alpha.f = 0.15), border = NA)
    graphics::polygon(c(fc_df$time, rev(fc_df$time)),
                      c(fc_df$upper_95, rev(fc_df$lower_95)),
                      col = grDevices::adjustcolor("steelblue", alpha.f = 0.25), border = NA)
    graphics::polygon(c(fc_df$time, rev(fc_df$time)),
                      c(fc_df$upper_90, rev(fc_df$lower_90)),
                      col = grDevices::adjustcolor("steelblue", alpha.f = 0.35), border = NA)

    # Forecast line
    graphics::lines(fc_df$time, fc_df$survival, col = "navy", lwd = lwd)
    graphics::lines(fc_df$time, fc_df$upper_95, col = "blue", lty = 3)
    graphics::lines(fc_df$time, fc_df$lower_95, col = "blue", lty = 3)

    graphics::legend("topright",
                     legend = c("Forecast Mean", "90% Interval", "95% Interval", "99% Interval"),
                     col = c("navy", grDevices::adjustcolor("steelblue", alpha.f = 0.5),
                             grDevices::adjustcolor("steelblue", alpha.f = 0.4),
                             grDevices::adjustcolor("steelblue", alpha.f = 0.2)),
                     lwd = c(lwd, 6, 6, 6), bg = "white", bty = "o")
  } else if (which == "qq") {
    r_quant <- x$residuals
    n_r <- length(r_quant)

    if (is.null(xlab)) xlab <- "Theoretical N(0,1) Quantiles"
    if (is.null(ylab)) ylab <- "Randomized Quantile Residuals"
    if (is.null(main)) main <- "Normal Q-Q Plot of Randomized Quantile Residuals"

    qq_res <- stats::qqnorm(r_quant, plot.it = FALSE)
    ord_x <- sort(qq_res$x)
    ord_y <- sort(qq_res$y)

    graphics::plot(ord_x, ord_y, pch = 19, col = grDevices::adjustcolor("navy", alpha.f = 0.7),
                   xlab = xlab, ylab = ylab, main = main, ...)
    graphics::grid(col = "gray90", lty = 1)
    graphics::abline(0, 1, col = "firebrick", lwd = 2)

    # 95% theoretical confidence envelope
    p_grid <- stats::ppoints(n_r)
    se_env <- sqrt(p_grid * (1 - p_grid) / n_r) / stats::dnorm(stats::qnorm(p_grid))
    upper_env <- stats::qnorm(p_grid) + 1.96 * se_env
    lower_env <- stats::qnorm(p_grid) - 1.96 * se_env

    graphics::lines(stats::qnorm(p_grid), upper_env, col = "gray50", lty = 2)
    graphics::lines(stats::qnorm(p_grid), lower_env, col = "gray50", lty = 2)
  } else if (which == "residuals") {
    r_quant <- x$residuals
    if (is.null(xlab)) xlab <- "Observation Index"
    if (is.null(ylab)) ylab <- "Randomized Quantile Residuals"
    if (is.null(main)) main <- "Residuals vs Index"

    graphics::plot(seq_along(r_quant), r_quant, pch = 19,
                   col = grDevices::adjustcolor("navy", alpha.f = 0.7),
                   xlab = xlab, ylab = ylab, main = main, ...)
    graphics::grid(col = "gray90", lty = 1)
    graphics::abline(h = 0, col = "firebrick", lwd = 1.5)
    graphics::abline(h = c(-2, 2), col = "gray40", lty = 2)
    graphics::abline(h = c(-3, 3), col = "red", lty = 3)
  } else if (which == "hazard") {
    h_vals <- predict(x, newdata = newdata, type = "hazard", times = times)
    h_curve <- if (is.matrix(h_vals)) colMeans(h_vals) else h_vals

    if (is.null(col)) col <- "darkorange"
    if (is.null(xlab)) xlab <- "Time"
    if (is.null(ylab)) ylab <- "Baseline Hazard Rate h0(t)"
    if (is.null(main)) main <- "Fitted Hazard Rate Function"

    graphics::plot(times, h_curve, type = "l", col = col, lwd = lwd,
                   xlab = xlab, ylab = ylab, main = main, ...)
    graphics::grid(col = "gray90", lty = 1)
    graphics::lines(times, h_curve, col = col, lwd = lwd)
  } else if (which == "density") {
    d_vals <- predict(x, newdata = newdata, type = "density", times = times)
    d_curve <- if (is.matrix(d_vals)) colMeans(d_vals) else d_vals

    if (is.null(col)) col <- "darkgreen"
    if (is.null(xlab)) xlab <- "Time"
    if (is.null(ylab)) ylab <- "Population Density fp(t)"
    if (is.null(main)) main <- "Fitted Population Density Function"

    graphics::plot(times, d_curve, type = "l", col = col, lwd = lwd,
                   xlab = xlab, ylab = ylab, main = main, ...)
    graphics::grid(col = "gray90", lty = 1)
    graphics::lines(times, d_curve, col = col, lwd = lwd)
  } else if (which == "activation_compare") {
    plot_activation_compare(x, times = times, col = col, lwd = lwd, ...)
  }

  invisible(NULL)
}

#' @rdname plot.ziprom_fit
#' @export
plot_survival_fit <- function(x, ...) {
  plot.ziprom_fit(x, which = "survival", ...)
}

#' @rdname plot.ziprom_fit
#' @export
plot_forecast <- function(x, ...) {
  plot.ziprom_fit(x, which = "forecast", ...)
}

#' @rdname plot.ziprom_fit
#' @export
plot_quantile_residuals <- function(x, ...) {
  plot.ziprom_fit(x, which = "qq", ...)
}

#' Compare Survival Curves Across Activation Schemes
#'
#' Fits and compares the fitted population survival curves under Minimum,
#' Maximum, and Random activation schemes for a given dataset.
#'
#' @param object An object of class \code{"ziprom_fit"}.
#' @param times Optional numeric vector of evaluation times.
#' @param col Colors for the three activation schemes.
#' @param lwd Line width (default is 2).
#' @param ... Additional graphical parameters.
#'
#' @return Invisible list of fitted models for each activation scheme.
#' @importFrom stats predict
#' @export
#' @examples
#' set.seed(42)
#' sim_dat <- sim_ziprom(n = 60, activation = "random")
#' fit <- ziprom_fit(survival::Surv(time, status) ~ x, data = sim_dat, activation = "random")
#' plot_activation_compare(fit)
plot_activation_compare <- function(object, times = NULL, col = NULL, lwd = 2, ...) {
  schemes <- c("random", "minimum", "maximum")
  if (is.null(col)) {
    col <- c("navy", "firebrick", "darkgreen")
  }

  max_t <- if (length(object$time) > 0) max(object$time) else 10
  if (is.null(times)) {
    times <- seq(0, max_t, length.out = 150)
  }

  fits <- list()
  s_curves <- list()

  for (s in schemes) {
    if (object$activation == s) {
      fits[[s]] <- object
    } else {
      # Refit under scheme s
      fits[[s]] <- tryCatch({
        ziprom_fit(
          time = object$time,
          status = object$status,
          data = object$data,
          formula_p0 = ~ 1,
          formula_p1 = ~ 1,
          formula_shape = ~ 1,
          formula_scale = ~ 1,
          activation = s,
          baseline = object$baseline,
          frailty = object$frailty
        )
      }, error = function(e) object)
    }
    s_curves[[s]] <- colMeans(predict(fits[[s]], type = "survival", times = times))
  }

  graphics::plot(times, s_curves[["random"]], type = "l", col = col[1], lwd = lwd,
                 ylim = c(0, 1), xlab = "Time", ylab = "Population Survival Probability",
                 main = "Activation Scheme Comparison", ...)
  graphics::grid(col = "gray90", lty = 1)
  graphics::lines(times, s_curves[["minimum"]], col = col[2], lwd = lwd, lty = 2)
  graphics::lines(times, s_curves[["maximum"]], col = col[3], lwd = lwd, lty = 3)

  if (length(object$time) > 0) {
    km_fit <- survival::survfit(survival::Surv(object$time, object$status) ~ 1)
    graphics::lines(km_fit, col = "gray40", lty = 4, lwd = 1.5, mark.time = FALSE)
    leg_txt <- c(paste0("Random (AIC=", round(fits[["random"]]$aic, 1), ")"),
                 paste0("Minimum (AIC=", round(fits[["minimum"]]$aic, 1), ")"),
                 paste0("Maximum (AIC=", round(fits[["maximum"]]$aic, 1), ")"),
                 "Kaplan-Meier")
    leg_col <- c(col, "gray40")
    leg_lty <- c(1, 2, 3, 4)
  } else {
    leg_txt <- c("Random", "Minimum", "Maximum")
    leg_col <- col
    leg_lty <- c(1, 2, 3)
  }

  graphics::legend("topright", legend = leg_txt, col = leg_col,
                   lty = leg_lty, lwd = c(rep(lwd, 3), 1.5), bg = "white", bty = "o")

  invisible(fits)
}
