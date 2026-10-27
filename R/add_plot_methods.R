#' Diagnostic and Visualization Plots for Additive Shared Frailty Models
#'
#' Generates MCMC trace plots, autocorrelation (ACF) plots, posterior density plots,
#' baseline hazard and survival curves, and goodness-of-fit diagnostic plots.
#'
#' @param x An object of class \code{"add_frailty"} or \code{"add_frailty_bayes"}.
#' @param type Character string specifying plot type: \code{"trace"}, \code{"acf"},
#'   \code{"density"}, \code{"survival"}, \code{"hazard"}, or \code{"gof"}.
#' @param params Optional character vector of parameter names to plot. If \code{NULL},
#'   plots the first 4 parameters.
#' @param max_time Maximum time horizon for survival/hazard plots (default 100).
#' @param ... Additional graphical parameters passed to plotting functions.
#'
#' @return Invisible \code{NULL}. Produces base R graphics.
#'
#' @export
#' @references
#' Pandey, A., Hanagal, D. D., and Tyagi, S. (2022). Generalised Lindley shared
#' additive frailty regression model for bivariate survival data.
#' \emph{Statistics in Transition new series}, 23(4), 161-176.
#' \doi{10.2478/stattrans-2022-0048}
#'
#' Tyagi, S., Pandey, A., Hanagal, D. D., and Chesneau, C. (2024). A study on
#' comparisons of additive regression frailty models to counter heterogeneity:
#' Bayesian strategies and case study. \emph{Communications in Statistics -
#' Simulation and Computation}. \doi{10.1080/03610918.2024.2391874}
#'
#' @examples
#' \donttest{
#' data(kidney_bivariate)
#' fit_mle <- fit_add_frailty_mle(
#'   data = kidney_bivariate,
#'   frailty_type = "wlindley", baseline_dist = "gw"
#' )
#' plot(fit_mle, type = "survival")
#' plot(fit_mle, type = "hazard")
#' }
plot.add_frailty <- function(x,
                             type = c("trace", "acf", "density", "survival", "hazard", "gof"),
                             params = NULL,
                             max_time = 100,
                             ...) {
  type <- match.arg(type)
  is_bayes <- inherits(x, "add_frailty_bayes")

  pars <- x$coefficients
  par_names <- names(pars)
  if (is.null(par_names)) par_names <- paste0("p", seq_along(pars))

  if (is.null(params)) {
    params <- par_names[1:min(4, length(par_names))]
  } else {
    params <- intersect(params, par_names)
    if (length(params) == 0) params <- par_names[1:min(4, length(par_names))]
  }

  n_plots <- length(params)

  if (type %in% c("trace", "acf", "density")) {
    if (!is_bayes) {
      stop(sprintf("Plot type '%s' is only available for Bayesian MCMC fitted objects.", type))
    }

    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(oldpar))

    n_rows <- ceiling(sqrt(n_plots))
    n_cols <- ceiling(n_plots / n_rows)
    graphics::par(mfrow = c(n_rows, n_cols), mar = c(3.5, 3.5, 2.5, 1), mgp = c(2, 0.7, 0))

    if (type == "trace") {
      for (p in params) {
        chain1 <- x$chains[[1]][, p]
        graphics::plot(chain1, type = "l", col = "blue",
                       xlab = "Iteration", ylab = "Value",
                       main = paste("Trace:", p), ...)
        if (length(x$chains) > 1) {
          chain2 <- x$chains[[2]][, p]
          graphics::lines(chain2, col = "darkorange", lty = 2)
          graphics::legend("topright", legend = c("Chain 1", "Chain 2"),
                           col = c("blue", "darkorange"), lty = c(1, 2), bty = "n", cex = 0.8)
        }
      }
    } else if (type == "acf") {
      for (p in params) {
        chain_vals <- x$posterior_samples[, p]
        stats::acf(chain_vals, main = paste("ACF:", p), col = "darkblue", ...)
      }
    } else if (type == "density") {
      for (p in params) {
        chain_vals <- x$posterior_samples[, p]
        dens <- stats::density(chain_vals)
        graphics::plot(dens, main = paste("Posterior Density:", p),
                       xlab = p, ylab = "Density", col = "darkgreen", lwd = 2, ...)
        graphics::polygon(dens, col = grDevices::rgb(0, 0.5, 0, 0.2), border = NA)
        graphics::abline(v = mean(chain_vals), col = "red", lty = 2, lwd = 2)
      }
    }
  } else if (type %in% c("survival", "hazard")) {
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(oldpar))

    graphics::par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
    t_grid <- seq(0.1, max_time, length.out = 200)

    z1 <- pars[1]; d1 <- pars[2]; x1 <- pars[3]
    z2 <- pars[4]; d2 <- pars[5]; x2 <- pars[6]

    if (x$model_info$baseline_dist == "gw") {
      s1 <- s_gw(t_grid, d1, z1, x1)
      s2 <- s_gw(t_grid, d2, z2, x2)
      h1 <- h_gw(t_grid, d1, z1, x1)
      h2 <- h_gw(t_grid, d2, z2, x2)
    } else {
      s1 <- s_gll2(t_grid, d1, z1, x1)
      s2 <- s_gll2(t_grid, d2, z2, x2)
      h1 <- h_gll2(t_grid, d1, z1, x1)
      h2 <- h_gll2(t_grid, d2, z2, x2)
    }

    if (type == "survival") {
      graphics::plot(t_grid, s1, type = "l", col = "navy", lwd = 2,
                     xlab = "Time (days)", ylab = "Survival S0(t)",
                     main = "Component 1 Baseline Survival", ylim = c(0, 1), ...)
      graphics::grid()
      graphics::plot(t_grid, s2, type = "l", col = "darkred", lwd = 2,
                     xlab = "Time (days)", ylab = "Survival S0(t)",
                     main = "Component 2 Baseline Survival", ylim = c(0, 1), ...)
      graphics::grid()
    } else {
      graphics::plot(t_grid, h1, type = "l", col = "navy", lwd = 2,
                     xlab = "Time (days)", ylab = "Hazard Rate h0(t)",
                     main = "Component 1 Baseline Hazard", ...)
      graphics::grid()
      graphics::plot(t_grid, h2, type = "l", col = "darkred", lwd = 2,
                     xlab = "Time (days)", ylab = "Hazard Rate h0(t)",
                     main = "Component 2 Baseline Hazard", ...)
      graphics::grid()
    }
  } else if (type == "gof") {
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(oldpar))

    graphics::par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
    t_grid <- seq(0.1, max_time, length.out = 200)

    z1 <- pars[1]; d1 <- pars[2]; x1 <- pars[3]
    z2 <- pars[4]; d2 <- pars[5]; x2 <- pars[6]

    if (x$model_info$baseline_dist == "gw") {
      f1 <- 1 - s_gw(t_grid, d1, z1, x1)
      f2 <- 1 - s_gw(t_grid, d2, z2, x2)
    } else {
      f1 <- 1 - s_gll2(t_grid, d1, z1, x1)
      f2 <- 1 - s_gll2(t_grid, d2, z2, x2)
    }

    graphics::plot(t_grid, f1, type = "l", col = "blue", lwd = 2,
                   xlab = "Time", ylab = "Cumulative Probability F(t)",
                   main = "GOF: Component 1 CDF", ylim = c(0, 1), ...)
    graphics::grid()
    graphics::plot(t_grid, f2, type = "l", col = "red", lwd = 2,
                   xlab = "Time", ylab = "Cumulative Probability F(t)",
                   main = "GOF: Component 2 CDF", ylim = c(0, 1), ...)
    graphics::grid()
  }

  invisible(NULL)
}
