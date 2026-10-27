#' Evaluation Metrics and Residual Diagnostics for Zero-Inflated Promotion Models
#'
#' Computes statistical goodness-of-fit and forecast accuracy measures including
#' Mean Squared Error (MSE), Root Mean Squared Error (RMSE), Mean Absolute Error (MAE),
#' Brier Score, and Harrell's Concordance Index (C-index), as well as exact
#' normalized randomized quantile residuals.
#'
#' @param object An object of class \code{"ziprom_fit"}.
#' @param newdata Optional new data frame to evaluate model performance on out-of-sample data.
#' @param times Optional numeric evaluation times for Brier Score / MSE.
#' @param ... Additional arguments.
#'
#' @return
#' For \code{evaluate_ziprom}: A list containing:
#' \describe{
#'   \item{\code{MSE}}{Mean Squared Error between observed status and fitted survival probabilities.}
#'   \item{\code{RMSE}}{Root Mean Squared Error (\eqn{\sqrt{\text{MSE}}}).}
#'   \item{\code{MAE}}{Mean Absolute Error.}
#'   \item{\code{Brier_Score}}{Integrated / average Brier score across time horizon.}
#'   \item{\code{C_index}}{Harrell's Concordance Index (discrimination metric).}
#'   \item{\code{LogLik}}{Maximum Log-Likelihood value.}
#'   \item{\code{AIC}}{Akaike Information Criterion.}
#'   \item{\code{BIC}}{Bayesian Information Criterion.}
#'   \item{\code{CAIC}}{Consistent Akaike Information Criterion.}
#'   \item{\code{HQIC}}{Hannan-Quinn Information Criterion.}
#' }
#'
#' For \code{ziprom_residuals}: A numeric vector of normalized randomized quantile residuals.
#'
#' @references
#' Dunn, P. K., and Smyth, G. K. (1996).
#' Randomized Quantile Residuals.
#' \emph{Journal of Computational and Graphical Statistics}, 5(3), 236-244. \doi{10.1080/10618600.1996.10474708}
#'
#' Assuncao, D., Ramos, P. L., Tomazella, V., Leao, J., Osses, L., and Rodrigues, A. (2026).
#' Zero-Inflated Promotion Cure-Rate Models With Frailty and Flexible Activation Schemes.
#' \emph{Biometrical Journal}, 68:e70162. \doi{10.1002/bimj.70162}
#'
#' @importFrom stats predict
#' @export
#' @examples
#' set.seed(42)
#' sim_dat <- sim_ziprom(n = 80, activation = "random")
#' fit <- ziprom_fit(survival::Surv(time, status) ~ x, data = sim_dat, activation = "random")
#' evals <- evaluate_ziprom(fit)
#' print(evals)
evaluate_ziprom <- function(object, newdata = NULL, times = NULL, ...) {
  if (is.null(newdata)) {
    time_vec <- object$time
    status_vec <- object$status
    b1 <- object$coef_list$beta1
    b2 <- object$coef_list$beta2
    b3 <- object$coef_list$beta3
    b4 <- object$coef_list$beta4
    alpha_val <- if (object$frailty) object$coef_list$alpha else 0

    probs <- compute_p0_p1(object$X1, object$X2, b1, b2)
    shape_vec <- exp(as.vector(object$X3 %*% b3))
    scale_vec <- exp(as.vector(object$X4 %*% b4))

    s_star <- ifelse(time_vec == 0, 1,
                     surv_proper(time_vec, probs$theta, alpha_val, shape_vec, scale_vec,
                                 object$activation, object$baseline))
    fitted_s_obs <- ifelse(time_vec == 0, 1 - probs$p0, probs$p1 + probs$susceptible * s_star)

    # Compute predicted survival at median follow-up time for discrimination
    t_med <- stats::median(time_vec[time_vec > 0], na.rm = TRUE)
    if (is.na(t_med) || t_med <= 0) t_med <- 1
    s_star_med <- surv_proper(t_med, probs$theta, alpha_val, shape_vec, scale_vec,
                              object$activation, object$baseline)
    s_med <- probs$p1 + probs$susceptible * s_star_med
  } else {
    newdata <- as.data.frame(newdata)
    if ("time" %in% names(newdata) && "status" %in% names(newdata)) {
      time_vec <- as.numeric(newdata$time)
      status_vec <- as.numeric(newdata$status)
    } else if (!is.null(object$call$formula)) {
      mf <- stats::model.frame(object$call$formula, data = newdata)
      y <- stats::model.response(mf)
      time_vec <- as.numeric(y[, 1])
      status_vec <- as.numeric(y[, 2])
    } else {
      time_vec <- object$time
      status_vec <- object$status
    }
    pred_s_mat <- predict(object, newdata = newdata, type = "survival", times = time_vec)
    fitted_s_obs <- if (is.matrix(pred_s_mat)) diag(pred_s_mat) else pred_s_mat

    t_med <- stats::median(time_vec[time_vec > 0], na.rm = TRUE)
    if (is.na(t_med) || t_med <= 0) t_med <- 1
    s_med_mat <- predict(object, newdata = newdata, type = "survival", times = t_med)
    s_med <- if (is.matrix(s_med_mat)) as.numeric(s_med_mat) else as.numeric(s_med_mat)
  }

  fitted_s_obs <- pmin(pmax(as.numeric(fitted_s_obs), 1e-12), 1 - 1e-12)
  n <- length(time_vec)

  # 1. MSE and MAE
  target_event <- 1 - status_vec
  mse_val <- mean((target_event - fitted_s_obs)^2, na.rm = TRUE)
  rmse_val <- sqrt(mse_val)
  mae_val <- mean(abs(target_event - fitted_s_obs), na.rm = TRUE)

  # 2. Time-dependent Brier Score across evaluation times
  if (is.null(times)) {
    pos_times <- time_vec[time_vec > 0 & status_vec == 1]
    if (length(pos_times) >= 5) {
      times_eval <- stats::quantile(pos_times, probs = seq(0.1, 0.9, length.out = 10), na.rm = TRUE)
    } else {
      times_eval <- seq(0.1, max(time_vec, na.rm = TRUE), length.out = 10)
    }
  } else {
    times_eval <- times
  }

  brier_scores <- numeric(length(times_eval))
  for (j in seq_along(times_eval)) {
    tj <- times_eval[j]
    if (is.null(newdata)) {
      s_star_j <- surv_proper(tj, probs$theta, alpha_val, shape_vec, scale_vec,
                              object$activation, object$baseline)
      s_tj <- probs$p1 + probs$susceptible * s_star_j
    } else {
      s_tj <- predict(object, newdata = newdata, type = "survival", times = tj)
      if (is.matrix(s_tj)) s_tj <- as.numeric(s_tj)
    }
    y_tj <- ifelse(time_vec > tj, 1, 0)
    brier_scores[j] <- mean((y_tj - s_tj)^2, na.rm = TRUE)
  }
  avg_brier <- mean(brier_scores, na.rm = TRUE)

  # 3. Harrell's Concordance Index (C-index)
  # Uses survival predictor (higher s_med implies longer expected survival)
  c_res <- tryCatch({
    cc <- survival::concordance(survival::Surv(time_vec, status_vec) ~ s_med)
    as.numeric(cc$concordance)
  }, error = function(e) {
    valid_idx <- !is.na(time_vec) & !is.na(status_vec) & !is.na(s_med)
    t_v <- time_vec[valid_idx]
    d_v <- status_vec[valid_idx]
    p_v <- s_med[valid_idx]
    n_v <- length(t_v)
    concordant <- 0
    discordant <- 0
    tied <- 0
    pairs <- 0
    if (n_v >= 2) {
      for (i in 1:(n_v - 1)) {
        for (k in (i + 1):n_v) {
          if (t_v[i] < t_v[k] && d_v[i] == 1) {
            pairs <- pairs + 1
            if (p_v[i] < p_v[k]) concordant <- concordant + 1
            else if (p_v[i] > p_v[k]) discordant <- discordant + 1
            else tied <- tied + 1
          } else if (t_v[k] < t_v[i] && d_v[k] == 1) {
            pairs <- pairs + 1
            if (p_v[k] < p_v[i]) concordant <- concordant + 1
            else if (p_v[k] > p_v[i]) discordant <- discordant + 1
            else tied <- tied + 1
          }
        }
      }
    }
    if (pairs > 0) (concordant + 0.5 * tied) / pairs else 0.5
  })

  list(
    MSE = round(mse_val, 6),
    RMSE = round(rmse_val, 6),
    MAE = round(mae_val, 6),
    Brier_Score = round(avg_brier, 6),
    C_index = round(c_res, 4),
    LogLik = round(object$loglik, 4),
    AIC = round(object$aic, 4),
    BIC = round(object$bic, 4),
    CAIC = round(object$caic, 4),
    HQIC = round(object$hqic, 4)
  )
}

#' @rdname evaluate_ziprom
#' @export
ziprom_residuals <- function(object, ...) {
  residuals.ziprom_fit(object, type = "quantile", ...)
}
