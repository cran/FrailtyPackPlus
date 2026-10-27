#' Prediction and Multi-Level Significance Forecasting for ziprom_fit Models
#'
#' Generates predicted survival probabilities, susceptible survival, cumulative hazard,
#' hazard rates, probability density, zero/cure proportions, and future forecasts
#' with 90\%, 95\%, and 99\% confidence / significance intervals.
#'
#' @param object An object of class \code{"ziprom_fit"}.
#' @param newdata Optional data frame containing covariate values for prediction.
#'   If omitted, predictions are made at the sample mean of the covariates.
#' @param type Character string indicating the type of prediction:
#'   \itemize{
#'     \item \code{"survival"} (default): Population survival probability \eqn{S_p(t)}.
#'     \item \code{"susceptible_survival"}: Proper survival for susceptible subjects \eqn{S_p^*(t)}.
#'     \item \code{"hazard"}: Baseline or population hazard rate.
#'     \item \code{"cumhazard"}: Cumulative baseline hazard \eqn{H_0(t)}.
#'     \item \code{"density"}: Population density \eqn{f_p(t)}.
#'     \item \code{"p0"}: Zero-failure probability \eqn{p_0}.
#'     \item \code{"p1"}: Cure fraction \eqn{p_1}.
#'     \item \code{"forecast"}: Detailed forecast table across time points with multi-level intervals.
#'     \item \code{"quantile"}: Survival quantiles (times) for susceptible individuals.
#'   }
#' @param times Numeric vector of future or evaluation time points.
#'   If \code{NULL}, a 100-point regular grid from 0 to \eqn{\max(t)} is constructed.
#' @param quantiles Numeric vector of quantile levels (in (0, 1)) for \code{type = "quantile"}.
#'   Default is \code{c(0.25, 0.50, 0.75)}.
#' @param level Numeric vector of confidence levels for interval estimation.
#'   Default is \code{c(0.90, 0.95, 0.99)}.
#' @param ... Additional arguments.
#'
#' @return
#' Depending on \code{type}:
#' \itemize{
#'   \item For \code{"survival"}, \code{"susceptible_survival"}, \code{"hazard"}, \code{"cumhazard"}, \code{"density"}:
#'     A numeric matrix (or vector) of predictions with rows corresponding to observations in \code{newdata}
#'     and columns to \code{times}.
#'   \item For \code{"p0"} or \code{"p1"}: A numeric vector of fitted probabilities.
#'   \item For \code{"quantile"}: A numeric matrix of estimated quantile survival times.
#'   \item For \code{"forecast"}: A comprehensive data frame containing \code{time}, \code{survival},
#'     \code{se}, \code{lower_90}, \code{upper_90}, \code{lower_95}, \code{upper_95},
#'     \code{lower_99}, and \code{upper_99}.
#' }
#'
#' @references
#' Assuncao, D., Ramos, P. L., Tomazella, V., Leao, J., Osses, L., and Rodrigues, A. (2026).
#' Zero-Inflated Promotion Cure-Rate Models With Frailty and Flexible Activation Schemes.
#' \emph{Biometrical Journal}, 68:e70162. \doi{10.1002/bimj.70162}
#'
#' @export
#' @examples
#' set.seed(42)
#' sim_dat <- sim_ziprom(n = 60, activation = "random")
#' fit <- ziprom_fit(survival::Surv(time, status) ~ x, data = sim_dat, activation = "random")
#' fc <- predict(fit, type = "forecast", times = seq(0, 5, length.out = 10))
#' print(head(fc))
predict.ziprom_fit <- function(object, newdata = NULL,
                               type = c("survival", "susceptible_survival", "hazard",
                                        "cumhazard", "density", "p0", "p1",
                                        "forecast", "quantile"),
                               times = NULL,
                               quantiles = c(0.25, 0.50, 0.75),
                               level = c(0.90, 0.95, 0.99), ...) {
  type <- match.arg(type)

  # Prepare times
  if (is.null(times)) {
    max_t <- if (length(object$time) > 0) max(object$time) else 10
    times <- seq(0, max_t, length.out = 100)
  }
  times <- sort(as.numeric(times))

  # Build Design Matrices for newdata
  if (is.null(newdata)) {
    X1_new <- matrix(colMeans(object$X1), nrow = 1, dimnames = list(NULL, colnames(object$X1)))
    X2_new <- matrix(colMeans(object$X2), nrow = 1, dimnames = list(NULL, colnames(object$X2)))
    X3_new <- matrix(colMeans(object$X3), nrow = 1, dimnames = list(NULL, colnames(object$X3)))
    X4_new <- matrix(colMeans(object$X4), nrow = 1, dimnames = list(NULL, colnames(object$X4)))
  } else {
    newdata <- as.data.frame(newdata)
    # Extract model terms if available
    try_matrix <- function(X_orig, nd) {
      cn <- colnames(X_orig)
      res <- matrix(0, nrow = nrow(nd), ncol = ncol(X_orig), dimnames = list(NULL, cn))
      res[, 1] <- 1 # Intercept
      for (j in seq_len(ncol(X_orig))) {
        col_name <- cn[j]
        if (col_name %in% names(nd)) {
          res[, j] <- as.numeric(nd[[col_name]])
        }
      }
      res
    }
    X1_new <- try_matrix(object$X1, newdata)
    X2_new <- try_matrix(object$X2, newdata)
    X3_new <- try_matrix(object$X3, newdata)
    X4_new <- try_matrix(object$X4, newdata)
  }

  n_new <- nrow(X1_new)
  b1 <- object$coef_list$beta1
  b2 <- object$coef_list$beta2
  b3 <- object$coef_list$beta3
  b4 <- object$coef_list$beta4
  alpha_val <- if (object$frailty) object$coef_list$alpha else 0

  probs <- compute_p0_p1(X1_new, X2_new, b1, b2)
  p0_vec <- probs$p0
  p1_vec <- probs$p1
  theta_vec <- probs$theta
  susc_vec <- probs$susceptible
  shape_vec <- exp(as.vector(X3_new %*% b3))
  scale_vec <- exp(as.vector(X4_new %*% b4))

  if (type == "p0") {
    return(p0_vec)
  }
  if (type == "p1") {
    return(p1_vec)
  }

  if (type == "quantile") {
    res_q <- matrix(NA_real_, nrow = n_new, ncol = length(quantiles),
                    dimnames = list(NULL, paste0("q_", quantiles * 100, "%")))
    for (i in 1:n_new) {
      for (k in seq_along(quantiles)) {
        res_q[i, k] <- quantile_proper(quantiles[k], theta_vec[i], alpha_val,
                                       shape_vec[i], scale_vec[i],
                                       object$activation, object$baseline)
      }
    }
    return(res_q)
  }

  if (type == "susceptible_survival") {
    res_mat <- matrix(NA_real_, nrow = n_new, ncol = length(times),
                      dimnames = list(NULL, paste0("t=", round(times, 3))))
    for (i in 1:n_new) {
      res_mat[i, ] <- ifelse(times == 0, 1,
                             surv_proper(times, theta_vec[i], alpha_val, shape_vec[i], scale_vec[i],
                                         object$activation, object$baseline))
    }
    return(res_mat)
  }

  if (type == "cumhazard") {
    res_mat <- matrix(NA_real_, nrow = n_new, ncol = length(times),
                      dimnames = list(NULL, paste0("t=", round(times, 3))))
    for (i in 1:n_new) {
      res_mat[i, ] <- base_cumhaz(times, shape_vec[i], scale_vec[i], object$baseline)
    }
    return(res_mat)
  }

  if (type == "hazard") {
    res_mat <- matrix(NA_real_, nrow = n_new, ncol = length(times),
                      dimnames = list(NULL, paste0("t=", round(times, 3))))
    for (i in 1:n_new) {
      res_mat[i, ] <- base_haz(times, shape_vec[i], scale_vec[i], object$baseline)
    }
    return(res_mat)
  }

  if (type == "density") {
    res_mat <- matrix(NA_real_, nrow = n_new, ncol = length(times),
                      dimnames = list(NULL, paste0("t=", round(times, 3))))
    for (i in 1:n_new) {
      res_mat[i, ] <- ziprom_core_dens(times, p0_vec[i], p1_vec[i], alpha_val,
                                       shape_vec[i], scale_vec[i],
                                       object$activation, object$baseline)
    }
    return(res_mat)
  }

  if (type %in% c("survival", "forecast")) {
    # Compute population survival matrix
    surv_mat <- matrix(NA_real_, nrow = n_new, ncol = length(times),
                       dimnames = list(NULL, paste0("t=", round(times, 3))))
    for (i in 1:n_new) {
      s_star_t <- ifelse(times == 0, 1,
                         surv_proper(times, theta_vec[i], alpha_val, shape_vec[i], scale_vec[i],
                                     object$activation, object$baseline))
      surv_mat[i, ] <- ifelse(times == 0, 1 - p0_vec[i], p1_vec[i] + susc_vec[i] * s_star_t)
    }

    if (type == "survival") {
      return(surv_mat)
    }

    # Detailed Forecasting with Multi-Level Significance Intervals
    # Compute pointwise delta method SE for first observation / representative profile
    x1_rep <- X1_new[1, , drop = FALSE]
    x2_rep <- X2_new[1, , drop = FALSE]
    x3_rep <- X3_new[1, , drop = FALSE]
    x4_rep <- X4_new[1, , drop = FALSE]

    vcov_mat <- object$vcov
    par_hat <- object$coefficients

    # Numerical gradient function for Sp(t)
    eval_surv_par <- function(p_vec, t_pt) {
      p1_l <- ncol(x1_rep)
      p2_l <- ncol(x2_rep)
      p3_l <- ncol(x3_rep)
      p4_l <- ncol(x4_rep)

      b1_p <- p_vec[1:p1_l]
      b2_p <- p_vec[(p1_l + 1):(p1_l + p2_l)]
      b3_p <- p_vec[(p1_l + p2_l + 1):(p1_l + p2_l + p3_l)]
      b4_p <- p_vec[(p1_l + p2_l + p3_l + 1):(p1_l + p2_l + p3_l + p4_l)]
      a_p <- if (object$frailty) exp(p_vec[length(p_vec)]) else 0

      pr <- compute_p0_p1(x1_rep, x2_rep, b1_p, b2_p)
      sh <- exp(as.vector(x3_rep %*% b3_p))
      sc <- exp(as.vector(x4_rep %*% b4_p))

      s_st <- if (t_pt == 0) 1 else surv_proper(t_pt, pr$theta, a_p, sh, sc, object$activation, object$baseline)
      if (t_pt == 0) 1 - pr$p0 else pr$p1 + pr$susceptible * s_st
    }

    se_forecast <- numeric(length(times))
    eps <- 1e-5
    for (k in seq_along(times)) {
      t_k <- times[k]
      grad_k <- numeric(length(par_hat))
      base_val <- eval_surv_par(par_hat, t_k)
      for (p in seq_along(par_hat)) {
        par_pert <- par_hat
        par_pert[p] <- par_pert[p] + eps
        grad_k[p] <- (eval_surv_par(par_pert, t_k) - base_val) / eps
      }
      se_forecast[k] <- sqrt(pmax(as.numeric(t(grad_k) %*% vcov_mat %*% grad_k), 1e-12))
    }

    s_fitted <- as.numeric(surv_mat[1, ])
    z90 <- stats::qnorm(0.95)
    z95 <- stats::qnorm(0.975)
    z99 <- stats::qnorm(0.995)

    forecast_df <- data.frame(
      time = times,
      survival = round(s_fitted, 6),
      se = round(se_forecast, 6),
      lower_90 = round(pmax(s_fitted - z90 * se_forecast, 0), 6),
      upper_90 = round(pmin(s_fitted + z90 * se_forecast, 1), 6),
      lower_95 = round(pmax(s_fitted - z95 * se_forecast, 0), 6),
      upper_95 = round(pmin(s_fitted + z95 * se_forecast, 1), 6),
      lower_99 = round(pmax(s_fitted - z99 * se_forecast, 0), 6),
      upper_99 = round(pmin(s_fitted + z99 * se_forecast, 1), 6)
    )
    return(forecast_df)
  }
}

#' Multi-Level Significance Survival Forecasting
#'
#' Helper function for extracting survival forecasts with 90\%, 95\%, and 99\%
#' level of significance intervals.
#'
#' @param object An object of class \code{"ziprom_fit"}.
#' @param newdata Optional data frame containing covariate values for forecasting.
#' @param times Optional numeric vector of forecast time points.
#' @param level Numeric vector of confidence levels. Default is \code{c(0.90, 0.95, 0.99)}.
#' @param ... Additional arguments passed to \code{\link{predict.ziprom_fit}}.
#'
#' @return A data frame of survival forecast probabilities with standard errors
#'   and multi-level significance intervals.
#' @export
#' @examples
#' set.seed(123)
#' dat <- sim_ziprom(n = 60, activation = "random")
#' fit <- ziprom_fit(survival::Surv(time, status) ~ x, data = dat, activation = "random")
#' fc_res <- forecast_ziprom(fit, times = c(0.5, 1.0, 2.0, 3.0, 5.0))
#' print(fc_res)
forecast_ziprom <- function(object, newdata = NULL, times = NULL,
                            level = c(0.90, 0.95, 0.99), ...) {
  predict.ziprom_fit(object, newdata = newdata, type = "forecast",
                     times = times, level = level, ...)
}
