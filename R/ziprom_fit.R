#' Fit Zero-Inflated Promotion Cure-Rate Models with Frailty
#'
#' Fits generalized zero-inflated promotion cure-rate regression models
#' with Gamma frailty and flexible activation schemes (minimum, maximum, or random)
#' using Maximum Likelihood Estimation (MLE).
#'
#' @param formula An object of class \code{\link[survival]{Surv}} or \code{\link[stats]{formula}},
#'   e.g. \code{survival::Surv(time, status) ~ x1 + x2}. If specified, the same linear predictor is used
#'   for components unless separate formulas are supplied.
#' @param data A data frame containing the variables named in the formulas.
#'   Missing values (\code{NA}) in the response or covariates are automatically removed
#'   prior to model fitting via listwise deletion (\code{\link[stats]{na.omit}}).
#' @param formula_p0 Optional one-sided formula for the zero-inflation component \eqn{p_0}, e.g. \code{~ x1}.
#' @param formula_p1 Optional one-sided formula for the cure fraction component \eqn{p_1}, e.g. \code{~ x2}.
#' @param formula_shape Optional one-sided formula for the baseline shape parameter \eqn{\alpha_1}, e.g. \code{~ 1}.
#' @param formula_scale Optional one-sided formula for the baseline scale parameter \eqn{\alpha_2}, e.g. \code{~ 1}.
#' @param time Optional numeric vector of observed survival times (used if \code{formula} is not provided).
#' @param status Optional numeric vector of event indicators (1 = event, 0 = censored).
#' @param activation Character string specifying the activation scheme:
#'   \code{"random"} (default), \code{"minimum"}, or \code{"maximum"}.
#' @param baseline Character string specifying the baseline activation-time distribution:
#'   \code{"weibull"} (default), \code{"gompertz"}, or \code{"exponential"}.
#' @param frailty Logical value (\code{TRUE} or \code{FALSE}) indicating whether Gamma frailty
#'   should be included in the model. Default is \code{TRUE}. If \code{FALSE}, fits the no-frailty model.
#' @param init Optional numeric vector of starting values for the optimization.
#' @param control Optional list of control parameters passed to \code{\link[stats]{optim}}.
#'
#' @details
#' Missing values (\code{NA}) in the response or covariates are omitted prior to estimation.
#' Observations with survival time equal to zero (\eqn{t = 0}) represent instantaneous failures and
#' are naturally accommodated via the zero-inflation probability mass \eqn{p_0}.
#'
#' @return An S3 object of class \code{"ziprom_fit"} containing:
#' \describe{
#'   \item{\code{coefficients}}{Vector of maximum likelihood estimates for all parameters.}
#'   \item{\code{coef_list}}{List of estimated coefficients partitioned by component (\code{beta1}, \code{beta2}, \code{beta3}, \code{beta4}, \code{alpha}).}
#'   \item{\code{se}}{Standard errors of parameter estimates obtained from the inverted observed Hessian.}
#'   \item{\code{z_values}}{Wald z-test statistics (\eqn{\hat{\beta} / \text{SE}}).}
#'   \item{\code{p_values}}{Two-sided asymptotic p-values for Wald tests.}
#'   \item{\code{ci_90}}{90\% Wald confidence intervals for parameters.}
#'   \item{\code{ci_95}}{95\% Wald confidence intervals for parameters.}
#'   \item{\code{ci_99}}{99\% Wald confidence intervals for parameters.}
#'   \item{\code{vcov}}{Estimated asymptotic variance-covariance matrix.}
#'   \item{\code{hessian}}{Observed information / Hessian matrix at the optimum.}
#'   \item{\code{loglik}}{Maximum log-likelihood value.}
#'   \item{\code{npar}}{Number of free estimated parameters.}
#'   \item{\code{aic}}{Akaike Information Criterion (\eqn{2k - 2\ell}).}
#'   \item{\code{bic}}{Bayesian Information Criterion (\eqn{k \log(n) - 2\ell}).}
#'   \item{\code{hqic}}{Hannan-Quinn Information Criterion.}
#'   \item{\code{caic}}{Consistent Akaike Information Criterion.}
#'   \item{\code{p0_mean}}{Estimated zero-failure proportion at mean covariates.}
#'   \item{\code{p0_mean_se}}{Delta-method standard error for mean \eqn{p_0}.}
#'   \item{\code{p1_mean}}{Estimated cure fraction at mean covariates.}
#'   \item{\code{p1_mean_se}}{Delta-method standard error for mean \eqn{p_1}.}
#'   \item{\code{p0_ci_95}}{95\% confidence interval for mean \eqn{p_0}.}
#'   \item{\code{p1_ci_95}}{95\% confidence interval for mean \eqn{p_1}.}
#'   \item{\code{fitted_p0}}{Fitted zero-inflation probabilities for each subject.}
#'   \item{\code{fitted_p1}}{Fitted cure probabilities for each subject.}
#'   \item{\code{residuals}}{Normalized randomized quantile residuals.}
#'   \item{\code{n}}{Total sample size.}
#'   \item{\code{n_zeros}}{Number of instantaneous zero failures (\eqn{t = 0}).}
#'   \item{\code{n_events}}{Number of positive observed events (\eqn{t > 0, \delta = 1}).}
#'   \item{\code{n_censored}}{Number of right-censored observations (\eqn{t > 0, \delta = 0}).}
#'   \item{\code{activation}}{Activation scheme used.}
#'   \item{\code{baseline}}{Baseline hazard distribution used.}
#'   \item{\code{frailty}}{Logical indicating whether frailty was modeled.}
#'   \item{\code{convergence}}{Convergence code returned by \code{optim} (0 indicates success).}
#'   \item{\code{data}}{Data frame containing the modeled variables.}
#'   \item{\code{call}}{The matched function call.}
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
#' sim_dat <- sim_ziprom(n = 80, beta1 = c(-1.5, 0.4), beta2 = c(-0.8, 0.6),
#'                       beta3 = c(0.2, 0.1), beta4 = c(1.0, 0.5), alpha = 0.25,
#'                       activation = "random", baseline = "weibull")
#' fit <- ziprom_fit(survival::Surv(time, status) ~ x, data = sim_dat,
#'                   activation = "random", frailty = TRUE)
#' summary(fit)
ziprom_fit <- function(formula = NULL, data = NULL,
                       formula_p0 = NULL, formula_p1 = NULL,
                       formula_shape = NULL, formula_scale = NULL,
                       time = NULL, status = NULL,
                       activation = c("random", "minimum", "maximum"),
                       baseline = c("weibull", "gompertz", "exponential"),
                       frailty = TRUE,
                       init = NULL,
                       control = list()) {
  call_matched <- match.call()
  activation <- match.arg(activation)
  baseline <- match.arg(baseline)
  frailty <- isTRUE(frailty)

  # 1. Parse Data and Formulas
  if (!is.null(formula)) {
    if (is.null(data)) {
      data <- environment(formula)
    }
    mf <- stats::model.frame(formula, data = data, na.action = stats::na.omit)
    y <- stats::model.response(mf)

    if (inherits(y, "Surv")) {
      time_vec <- as.numeric(y[, 1])
      status_vec <- as.numeric(y[, 2])
    } else {
      time_vec <- as.numeric(y)
      if (is.null(status)) {
        stop("When response in formula is not a Surv object, 'status' vector must be provided.")
      }
      status_vec <- as.numeric(status)
    }

    # Default submodel formulas inherit RHS of formula
    rhs_formula <- stats::as.formula(paste("~", paste(attr(stats::terms(formula), "term.labels"), collapse = " + ")))
    if (length(attr(stats::terms(formula), "term.labels")) == 0) {
      rhs_formula <- ~ 1
    }

    f_p0 <- if (!is.null(formula_p0)) formula_p0 else rhs_formula
    f_p1 <- if (!is.null(formula_p1)) formula_p1 else rhs_formula
    f_shape <- if (!is.null(formula_shape)) formula_shape else ~ 1
    f_scale <- if (!is.null(formula_scale)) formula_scale else rhs_formula

    mf_p0 <- stats::model.frame(f_p0, data = mf)
    mf_p1 <- stats::model.frame(f_p1, data = mf)
    mf_shape <- stats::model.frame(f_shape, data = mf)
    mf_scale <- stats::model.frame(f_scale, data = mf)

    X1 <- stats::model.matrix(f_p0, mf_p0)
    X2 <- stats::model.matrix(f_p1, mf_p1)
    X3 <- stats::model.matrix(f_shape, mf_shape)
    X4 <- stats::model.matrix(f_scale, mf_scale)

    df_used <- mf
  } else {
    if (is.null(time) || is.null(status)) {
      stop("Either 'formula' or both 'time' and 'status' must be provided.")
    }
    time_vec <- as.numeric(time)
    status_vec <- as.numeric(status)

    valid_idx <- !(is.na(time_vec) | is.na(status_vec))
    time_vec <- time_vec[valid_idx]
    status_vec <- status_vec[valid_idx]

    if (is.null(data)) {
      X1 <- matrix(1, nrow = length(time_vec), ncol = 1, dimnames = list(NULL, "(Intercept)"))
      X2 <- matrix(1, nrow = length(time_vec), ncol = 1, dimnames = list(NULL, "(Intercept)"))
      X3 <- matrix(1, nrow = length(time_vec), ncol = 1, dimnames = list(NULL, "(Intercept)"))
      X4 <- matrix(1, nrow = length(time_vec), ncol = 1, dimnames = list(NULL, "(Intercept)"))
      df_used <- data.frame(time = time_vec, status = status_vec)
    } else {
      df_used <- data[valid_idx, , drop = FALSE]
      f_p0 <- if (!is.null(formula_p0)) formula_p0 else ~ 1
      f_p1 <- if (!is.null(formula_p1)) formula_p1 else ~ 1
      f_shape <- if (!is.null(formula_shape)) formula_shape else ~ 1
      f_scale <- if (!is.null(formula_scale)) formula_scale else ~ 1

      X1 <- stats::model.matrix(f_p0, stats::model.frame(f_p0, df_used))
      X2 <- stats::model.matrix(f_p1, stats::model.frame(f_p1, df_used))
      X3 <- stats::model.matrix(f_shape, stats::model.frame(f_shape, df_used))
      X4 <- stats::model.matrix(f_scale, stats::model.frame(f_scale, df_used))
    }
  }

  n <- length(time_vec)
  n_zeros <- sum(time_vec == 0)
  n_events <- sum(time_vec > 0 & status_vec == 1)
  n_censored <- sum(time_vec > 0 & status_vec == 0)

  # Parameter names
  name_p0 <- paste0("p0_", colnames(X1))
  name_p1 <- paste0("p1_", colnames(X2))
  name_shape <- paste0("shape_", colnames(X3))
  name_scale <- paste0("scale_", colnames(X4))

  param_names <- c(name_p0, name_p1, name_shape, name_scale)
  if (frailty) {
    param_names <- c(param_names, "log_alpha")
  }

  # 2. Smart Initial Values
  if (is.null(init)) {
    prop_zeros <- max(n_zeros / n, 0.05)
    prop_cure <- max(n_censored / n, 0.10)
    prop_susc <- max(1 - prop_zeros - prop_cure, 0.10)

    beta1_init <- rep(0, ncol(X1))
    beta1_init[1] <- log(prop_zeros / prop_susc)

    beta2_init <- rep(0, ncol(X2))
    beta2_init[1] <- log(prop_cure / prop_susc)

    pos_times <- time_vec[time_vec > 0]
    med_time <- if (length(pos_times) > 0) stats::median(pos_times) else 1.0
    med_time <- max(med_time, 0.1)

    beta3_init <- rep(0, ncol(X3)) # shape = exp(0) = 1
    beta4_init <- rep(0, ncol(X4))
    beta4_init[1] <- log(med_time / (log(2))^(1 / 1)) # scale

    init_par <- c(beta1_init, beta2_init, beta3_init, beta4_init)
    if (frailty) {
      init_par <- c(init_par, log(0.25))
    }
  } else {
    init_par <- init
    if (length(init_par) != length(param_names)) {
      stop(sprintf("Length of 'init' (%d) does not match required parameters (%d).",
                   length(init_par), length(param_names)))
    }
  }

  # 3. Optimization
  ctrl <- list(maxit = 2000, reltol = 1e-10)
  if (length(control) > 0) {
    ctrl[names(control)] <- control
  }

  opt_res <- stats::optim(
    par = init_par,
    fn = ziprom_loglik,
    time = time_vec,
    status = status_vec,
    X1 = X1,
    X2 = X2,
    X3 = X3,
    X4 = X4,
    activation = activation,
    baseline = baseline,
    frailty = frailty,
    method = "BFGS",
    hessian = TRUE,
    control = ctrl
  )

  # Fallback to Nelder-Mead if BFGS failed to converge
  if (opt_res$convergence != 0) {
    opt_res_nm <- stats::optim(
      par = opt_res$par,
      fn = ziprom_loglik,
      time = time_vec,
      status = status_vec,
      X1 = X1,
      X2 = X2,
      X3 = X3,
      X4 = X4,
      activation = activation,
      baseline = baseline,
      frailty = frailty,
      method = "Nelder-Mead",
      hessian = TRUE,
      control = ctrl
    )
    if (opt_res_nm$value < opt_res$value) {
      opt_res <- opt_res_nm
    }
  }

  par_hat <- opt_res$par
  names(par_hat) <- param_names
  loglik_val <- -opt_res$value
  k_par <- length(par_hat)

  # 4. Variance-Covariance Matrix and Standard Errors
  H <- opt_res$hessian
  dimnames(H) <- list(param_names, param_names)

  # Regularized inversion to ensure non-negative variance
  vcov_mat <- tryCatch({
    solve(H)
  }, error = function(e) {
    eig <- eigen(H, symmetric = TRUE)
    eig$values <- pmax(eig$values, 1e-6)
    eig$vectors %*% diag(1 / eig$values) %*% t(eig$vectors)
  })
  dimnames(vcov_mat) <- list(param_names, param_names)

  se_vec <- sqrt(pmax(diag(vcov_mat), 1e-10))
  names(se_vec) <- param_names

  z_vals <- par_hat / se_vec
  p_vals <- 2 * (1 - stats::pnorm(abs(z_vals)))

  # Confidence Intervals for parameters
  z90 <- stats::qnorm(0.95)
  z95 <- stats::qnorm(0.975)
  z99 <- stats::qnorm(0.995)

  ci90 <- cbind(Lower = par_hat - z90 * se_vec, Upper = par_hat + z90 * se_vec)
  ci95 <- cbind(Lower = par_hat - z95 * se_vec, Upper = par_hat + z95 * se_vec)
  ci99 <- cbind(Lower = par_hat - z99 * se_vec, Upper = par_hat + z99 * se_vec)

  # 5. Parameter List & Delta-Method for Proportions
  p1_len <- ncol(X1)
  p2_len <- ncol(X2)
  p3_len <- ncol(X3)
  p4_len <- ncol(X4)

  b1_hat <- par_hat[1:p1_len]
  b2_hat <- par_hat[(p1_len + 1):(p1_len + p2_len)]
  b3_hat <- par_hat[(p1_len + p2_len + 1):(p1_len + p2_len + p3_len)]
  b4_hat <- par_hat[(p1_len + p2_len + p3_len + 1):(p1_len + p2_len + p3_len + p4_len)]

  alpha_hat <- if (frailty) exp(par_hat[k_par]) else 0
  alpha_se <- if (frailty) alpha_hat * se_vec[k_par] else 0

  coef_list <- list(
    beta1 = b1_hat,
    beta2 = b2_hat,
    beta3 = b3_hat,
    beta4 = b4_hat,
    alpha = alpha_hat,
    alpha_se = alpha_se
  )

  # Proportion estimates at sample mean covariates
  x1_bar <- colMeans(X1)
  x2_bar <- colMeans(X2)
  eta1_bar <- sum(x1_bar * b1_hat)
  eta2_bar <- sum(x2_bar * b2_hat)
  denom_bar <- 1 + exp(eta1_bar) + exp(eta2_bar)

  p0_mean <- exp(eta1_bar) / denom_bar
  p1_mean <- exp(eta2_bar) / denom_bar

  # Gradients for p0 and p1
  grad_p0 <- c(p0_mean * (1 - p0_mean) * x1_bar, -p0_mean * p1_mean * x2_bar)
  grad_p1 <- c(-p0_mean * p1_mean * x1_bar, p1_mean * (1 - p1_mean) * x2_bar)

  vcov_b12 <- vcov_mat[1:(p1_len + p2_len), 1:(p1_len + p2_len), drop = FALSE]
  p0_mean_se <- as.numeric(sqrt(pmax(t(grad_p0) %*% vcov_b12 %*% grad_p0, 1e-12)))
  p1_mean_se <- as.numeric(sqrt(pmax(t(grad_p1) %*% vcov_b12 %*% grad_p1, 1e-12)))

  p0_ci_95 <- c(max(0, p0_mean - z95 * p0_mean_se), min(1, p0_mean + z95 * p0_mean_se))
  p1_ci_95 <- c(max(0, p1_mean - z95 * p1_mean_se), min(1, p1_mean + z95 * p1_mean_se))

  # Fitted probabilities per individual
  fitted_probs <- compute_p0_p1(X1, X2, b1_hat, b2_hat)
  fitted_p0 <- fitted_probs$p0
  fitted_p1 <- fitted_probs$p1

  # 6. Information Criteria
  aic_val <- 2 * k_par - 2 * loglik_val
  bic_val <- k_par * log(n) - 2 * loglik_val
  caic_val <- k_par * (log(n) + 1) - 2 * loglik_val
  hqic_val <- 2 * k_par * log(log(n)) - 2 * loglik_val

  # 7. Normalized Randomized Quantile Residuals
  shape_vec <- exp(as.vector(X3 %*% b3_hat))
  scale_vec <- exp(as.vector(X4 %*% b4_hat))
  theta_vec <- fitted_probs$theta
  susc_vec <- fitted_probs$susceptible

  s_star_fitted <- ifelse(time_vec == 0, 1,
                          surv_proper(time_vec, theta_vec, alpha_hat, shape_vec, scale_vec, activation, baseline))
  surv_pop_fitted <- ifelse(time_vec == 0, 1 - fitted_p0, fitted_p1 + susc_vec * s_star_fitted)
  surv_pop_fitted <- pmin(pmax(surv_pop_fitted, 1e-12), 1 - 1e-12)

  # Exact randomized quantile residual computation
  u_rand <- numeric(n)
  for (i in 1:n) {
    if (time_vec[i] == 0) {
      u_rand[i] <- stats::runif(1, 0, fitted_p0[i])
    } else if (status_vec[i] == 1) {
      u_rand[i] <- 1 - surv_pop_fitted[i]
    } else {
      u_rand[i] <- stats::runif(1, 1 - surv_pop_fitted[i], 1)
    }
  }
  u_rand <- pmin(pmax(u_rand, 1e-10), 1 - 1e-10)
  res_quant <- stats::qnorm(u_rand)

  structure(
    list(
      coefficients = par_hat,
      coef_list = coef_list,
      se = se_vec,
      z_values = z_vals,
      p_values = p_vals,
      ci_90 = ci90,
      ci_95 = ci95,
      ci_99 = ci99,
      vcov = vcov_mat,
      hessian = H,
      loglik = loglik_val,
      npar = k_par,
      aic = aic_val,
      bic = bic_val,
      hqic = hqic_val,
      caic = caic_val,
      p0_mean = p0_mean,
      p0_mean_se = p0_mean_se,
      p1_mean = p1_mean,
      p1_mean_se = p1_mean_se,
      p0_ci_95 = p0_ci_95,
      p1_ci_95 = p1_ci_95,
      fitted_p0 = fitted_p0,
      fitted_p1 = fitted_p1,
      residuals = res_quant,
      time = time_vec,
      status = status_vec,
      X1 = X1,
      X2 = X2,
      X3 = X3,
      X4 = X4,
      n = n,
      n_zeros = n_zeros,
      n_events = n_events,
      n_censored = n_censored,
      activation = activation,
      baseline = baseline,
      frailty = frailty,
      convergence = opt_res$convergence,
      data = df_used,
      call = call_matched
    ),
    class = "ziprom_fit"
  )
}

#' @importFrom survival Surv
#' @export
survival::Surv
