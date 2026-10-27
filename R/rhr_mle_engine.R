# Maximum Likelihood Estimation Engine for Shared RHR Frailty Models

#' Maximum Likelihood Estimation for Shared RHR Frailty Models
#'
#' Fits shared frailty regression models for bivariate left-censored survival data
#' based on reversed hazard rates using maximum likelihood estimation.
#'
#' @details
#' Maximizes the four-part left-censored log-likelihood function using the BFGS algorithm
#' on unconstrained log-transformed parameter space to guarantee positive parameters.
#' Standard errors are computed from the numerical Hessian using the multivariate Delta Method.
#' Missing values (\code{NA} or \code{NaN}) in survival times or censoring indicators are not permitted.
#'
#' @param formula An optional object of class \code{"formula"} of the form \code{cbind(t1, t2, d1, d2) ~ x1 + x2}.
#' @param data An optional data frame containing variables.
#' @param t1 Numeric vector of first failure/censoring times. Must be positive.
#' @param t2 Numeric vector of second failure/censoring times. Must be positive.
#' @param d1 Numeric or integer vector of censoring indicators for first component (\code{1} = observed, \code{0} = left-censored).
#' @param d2 Numeric or integer vector of censoring indicators for second component (\code{1} = observed, \code{0} = left-censored).
#' @param covars Optional numeric matrix of covariates.
#' @param frailty Character string specifying frailty distribution. Default is \code{"WL"}.
#' @param baseline Character string specifying baseline distribution. Default is \code{"GE"}.
#' @param init Optional numeric vector of initial parameters.
#' @param hessian Logical indicating whether to compute observed Fisher information Hessian. Default is \code{TRUE}.
#' @param control Optional list of control parameters passed to \code{\link[stats]{optim}}.
#'
#' @return An object of class \code{"rhr_mle_fit"} with the following components:
#' \describe{
#'   \item{\code{summary}}{A data frame containing parameter estimates, standard errors, z-scores, p-values, and Wald 95\% confidence intervals.}
#'   \item{\code{coefficients}}{A named numeric vector of parameter estimates.}
#'   \item{\code{se}}{A named numeric vector of asymptotic standard errors computed via the Delta Method.}
#'   \item{\code{vcov}}{The estimated asymptotic variance-covariance matrix.}
#'   \item{\code{loglik}}{The maximized log-likelihood value.}
#'   \item{\code{criteria}}{A list containing AIC, BIC, AICc, and HQIC information criteria.}
#'   \item{\code{kendall_tau}}{Numeric scalar representing estimated Kendall's tau.}
#'   \item{\code{frailty_variance}}{Numeric scalar representing estimated variance of the frailty variable.}
#'   \item{\code{opt}}{The raw optimization result returned by \code{\link[stats]{optim}}.}
#'   \item{\code{model_info}}{A list containing model specifications: frailty, baseline, sample size, number of covariates, covariate names, and parameter names.}
#'   \item{\code{data}}{A list containing the input data vectors and covariate matrix.}
#'   \item{\code{call}}{The matched call.}
#' }
#' @export
#'
#' @examples
#' data(austrian_twin)
#' fit <- fit_rhr_frailty_mle(
#'   cbind(t1, t2, d1, d2) ~ age + k2 + k3,
#'   data = austrian_twin[1:50, ],
#'   frailty = "WL", baseline = "GE"
#' )
#' print(fit)
#' summary(fit)
fit_rhr_frailty_mle <- function(formula = NULL, data = NULL, t1 = NULL, t2 = NULL, d1 = NULL, d2 = NULL,
                                covars = NULL, frailty = c("WL", "GL", "GIG", "Gamma", "IG", "Lindley", "None"),
                                baseline = c("GE", "EG", "GIR", "MIW", "GLL-I", "GLL-II", "GIE", "Weibull", "Exponential", "Log-Logistic"),
                                init = NULL, hessian = TRUE, control = list()) {
  frailty <- match.arg(frailty)
  baseline <- match.arg(baseline)
  
  if (!is.null(formula)) {
    if (is.null(data)) data <- environment(formula)
    mf <- stats::model.frame(formula, data)
    lhs <- stats::model.response(mf)
    t1 <- as.numeric(lhs[, 1])
    t2 <- as.numeric(lhs[, 2])
    d1 <- as.numeric(lhs[, 3])
    d2 <- as.numeric(lhs[, 4])
    covars <- stats::model.matrix(formula, mf)
    if (ncol(covars) > 1 && colnames(covars)[1] == "(Intercept)") {
      covars <- covars[, -1, drop = FALSE]
    }
  } else {
    t1 <- as.numeric(t1)
    t2 <- as.numeric(t2)
    d1 <- as.numeric(d1)
    d2 <- as.numeric(d2)
    if (!is.null(covars)) covars <- as.matrix(covars)
  }
  
  if (any(is.na(t1)) || any(is.na(t2)) || any(is.na(d1)) || any(is.na(d2))) {
    stop("Missing values (NA or NaN) detected in survival times or censoring indicators.")
  }
  
  n <- length(t1)
  p_cov <- if (!is.null(covars)) ncol(as.matrix(covars)) else 0
  cov_names <- if (p_cov > 0) colnames(covars) else character(0)
  if (p_cov > 0 && is.null(cov_names)) cov_names <- paste0("beta", seq_len(p_cov))
  
  b_names <- switch(baseline,
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
  
  p1_names <- paste0(b_names, "1")
  p2_names <- paste0(b_names, "2")
  
  f_names <- switch(frailty,
    "WL" = "theta",
    "GL" = c("theta", "mu"),
    "GIG" = c("theta", "theta2"),
    "Gamma" = "theta",
    "IG" = "theta",
    "Lindley" = "theta",
    "None" = character(0)
  )
  
  all_param_names <- c(p1_names, p2_names, f_names, cov_names)
  n_params <- length(all_param_names)
  n_pos <- n_params - p_cov
  
  if (is.null(init)) {
    init_p1 <- switch(baseline,
      "GE" = c(zeta = 10.0, delta = 0.10),
      "EG" = c(zeta = 8.0, delta = 0.10),
      "GIR" = c(zeta = 3.5, delta = 900.0),
      "MIW" = c(zeta = 3.0, delta = 1.0, xi = 0.05),
      "GLL-I" = c(lambda = 0.05, alpha = 30.0, gamma = 1.5),
      "GLL-II" = c(lambda = 0.05, alpha = 0.15, gamma = 1.5),
      "GIE" = c(zeta = 3.0, delta = 10.0),
      "Weibull" = c(zeta = 1.5, delta = 0.05),
      "Exponential" = c(delta = 0.05),
      "Log-Logistic" = c(lambda = 0.05, gamma = 1.5)
    )
    init_p2 <- init_p1
    init_f <- switch(frailty,
      "WL" = c(theta = 3.0),
      "GL" = c(theta = 3.0, mu = 1.5),
      "GIG" = c(theta = 2.0, theta2 = 2.0),
      "Gamma" = c(theta = 2.0),
      "IG" = c(theta = 2.0),
      "Lindley" = c(theta = 2.0),
      "None" = numeric(0)
    )
    init_beta <- if (p_cov > 0) rep(0.0, p_cov) else numeric(0)
    init <- c(init_p1, init_p2, init_f, init_beta)
    names(init) <- all_param_names
  }
  
  uncon_init <- init
  uncon_init[1:n_pos] <- log(pmax(init[1:n_pos], 1e-6))
  
  obj_fn <- function(uncon_p) {
    p <- uncon_p
    p[1:n_pos] <- exp(uncon_p[1:n_pos])
    names(p) <- all_param_names
    ll <- loglik_biv_rhr(p, t1, t2, d1, d2, covars, frailty, baseline)
    if (!is.finite(ll)) return(1e10)
    -ll
  }
  
  opt_res <- stats::optim(uncon_init, obj_fn, method = "BFGS", hessian = hessian, control = control)
  
  est_uncon <- opt_res$par
  est_params <- est_uncon
  est_params[1:n_pos] <- exp(est_uncon[1:n_pos])
  names(est_params) <- all_param_names
  
  max_ll <- -opt_res$value
  
  vcov_mat <- matrix(NA_real_, nrow = n_params, ncol = n_params, dimnames = list(all_param_names, all_param_names))
  se_vec <- rep(NA_real_, n_params)
  names(se_vec) <- all_param_names
  
  if (hessian && !inherits(try(solve(opt_res$hessian), silent = TRUE), "try-error")) {
    inv_hess <- solve(opt_res$hessian)
    diag_J <- c(est_params[1:n_pos], rep(1, p_cov))
    J <- diag(diag_J, nrow = n_params, ncol = n_params)
    vcov_mat <- J %*% inv_hess %*% J
    dimnames(vcov_mat) <- list(all_param_names, all_param_names)
    se_vec <- sqrt(pmax(diag(vcov_mat), 0))
  }
  
  z_vals <- est_params / pmax(se_vec, 1e-12)
  p_vals <- 2 * (1 - stats::pnorm(abs(z_vals)))
  ci_lower <- est_params - 1.96 * se_vec
  ci_upper <- est_params + 1.96 * se_vec
  ci_lower[1:n_pos] <- pmax(ci_lower[1:n_pos], 0)
  
  AIC_val <- 2 * n_params - 2 * max_ll
  BIC_val <- n_params * log(n) - 2 * max_ll
  AICc_val <- AIC_val + (2 * n_params * (n_params + 1)) / max(1, n - n_params - 1)
  HQIC_val <- 2 * n_params * log(log(n)) - 2 * max_ll
  
  th_est <- if (frailty != "None") est_params["theta"] else 1.0
  mu_est <- if (frailty == "GL") est_params["mu"] else NULL
  th2_est <- if (frailty == "GIG") est_params["theta2"] else NULL
  
  tau_est <- kendall_tau_rhr(frailty, th_est, mu_est, th2_est)
  var_est <- var_frailty(frailty, th_est, mu_est, th2_est)
  
  summary_table <- data.frame(
    Parameter = all_param_names,
    Estimate = est_params,
    StdError = se_vec,
    Z_value = z_vals,
    P_value = p_vals,
    CI_Lower = ci_lower,
    CI_Upper = ci_upper,
    stringsAsFactors = FALSE
  )
  
  out <- list(
    summary = summary_table,
    coefficients = est_params,
    se = se_vec,
    vcov = vcov_mat,
    loglik = max_ll,
    criteria = list(
      AIC = AIC_val,
      BIC = BIC_val,
      AICc = AICc_val,
      HQIC = HQIC_val
    ),
    kendall_tau = tau_est,
    frailty_variance = var_est,
    opt = opt_res,
    model_info = list(
      frailty = frailty,
      baseline = baseline,
      n_obs = n,
      p_cov = p_cov,
      cov_names = cov_names,
      param_names = all_param_names
    ),
    data = list(t1 = t1, t2 = t2, d1 = d1, d2 = d2, covars = covars),
    call = match.call()
  )
  
  class(out) <- "rhr_mle_fit"
  return(out)
}
