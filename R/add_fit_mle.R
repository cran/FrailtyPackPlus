#' Maximum Likelihood Estimation for Additive Shared Frailty Models
#'
#' Estimates unknown parameters of additive shared frailty models under bivariate
#' random censoring using numerical optimization of the log-likelihood with
#' transformation for bounded parameters and Hessian-based standard errors.
#'
#' @param formula Model formula for covariates (e.g. \code{~ age + sex + gn + an + pkd}).
#'   If \code{NULL}, an intercept-only or single-covariate model is fitted.
#' @param data A data frame containing survival times, censoring indicators, and covariates.
#' @param t1_col Character string naming the first survival time column in \code{data}.
#' @param delta1_col Character string naming the first censoring indicator column in \code{data}.
#' @param t2_col Character string naming the second survival time column in \code{data}.
#' @param delta2_col Character string naming the second censoring indicator column in \code{data}.
#' @param frailty_type Character string specifying the frailty distribution:
#'   \code{"genlindley"} (Generalized Lindley), \code{"wlindley"} (Weighted Lindley),
#'   \code{"gamma"} (Gamma), or \code{"invgaussian"} (Inverse Gaussian).
#' @param baseline_dist Character string specifying the baseline distribution:
#'   \code{"gw"} (Generalized Weibull) or \code{"gll2"} (Generalized Log-Logistic Type-II).
#' @param init Optional numeric vector of initial parameter values.
#' @param method Optimization method passed to \code{stats::optim} (default \code{"Nelder-Mead"}).
#' @param hessian Logical; if \code{TRUE}, computes the numerical Hessian for standard errors.
#' @param control List of control parameters passed to \code{stats::optim}.
#'
#' @return An object of class \code{c("add_frailty_mle", "add_frailty")} containing:
#' \describe{
#'   \item{coefficients}{Named vector of Maximum Likelihood parameter estimates.}
#'   \item{estimates}{Data frame containing estimates, standard errors, z-values, p-values, and 95\% confidence bounds.}
#'   \item{vcov}{Asymptotic covariance matrix of the estimates.}
#'   \item{loglik}{Maximized log-likelihood value.}
#'   \item{criteria}{Named list of AIC, BIC, AICc, HQIC values.}
#'   \item{kendall_tau}{Kendall's tau dependence measure evaluated at the MLE.}
#'   \item{frailty_variance}{Theoretical frailty variance evaluated at the MLE.}
#'   \item{convergence}{Convergence code from \code{stats::optim}.}
#'   \item{model_info}{List containing model specifications and call.}
#' }
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
#' sub_data <- kidney_bivariate[1:20, ]
#' fit_mle <- fit_add_frailty_mle(
#'   formula = ~ age,
#'   data = sub_data,
#'   frailty_type = "wlindley", baseline_dist = "gw",
#'   control = list(maxit = 100)
#' )
#' print(fit_mle)
#' }
fit_add_frailty_mle <- function(formula = NULL,
                                data,
                                t1_col = "t1",
                                delta1_col = "delta1",
                                t2_col = "t2",
                                delta2_col = "delta2",
                                frailty_type = c("genlindley", "wlindley", "gamma", "invgaussian"),
                                baseline_dist = c("gw", "gll2"),
                                init = NULL,
                                method = "Nelder-Mead",
                                hessian = TRUE,
                                control = list(maxit = 5000)) {
  call <- match.call()
  frailty_type <- match.arg(frailty_type)
  baseline_dist <- match.arg(baseline_dist)

  if (!is.data.frame(data)) {
    stop("'data' must be a data frame.")
  }

  t1 <- as.numeric(data[[t1_col]])
  delta1 <- as.numeric(data[[delta1_col]])
  t2 <- as.numeric(data[[t2_col]])
  delta2 <- as.numeric(data[[delta2_col]])

  if (is.null(formula)) {
    cov_cols <- setdiff(colnames(data), c(t1_col, delta1_col, t2_col, delta2_col, "id", "w"))
    if (length(cov_cols) > 0) {
      X <- as.matrix(data[, cov_cols, drop = FALSE])
    } else {
      X <- matrix(0, nrow = length(t1), ncol = 0)
    }
  } else {
    X <- stats::model.matrix(formula, data)
    if ("(Intercept)" %in% colnames(X)) {
      X <- X[, colnames(X) != "(Intercept)", drop = FALSE]
    }
  }

  n_cov <- ncol(X)
  cov_names <- if (n_cov > 0) colnames(X) else character(0)

  if (frailty_type == "genlindley") {
    par_names <- c("zeta1", "delta1", "xi1", "zeta2", "delta2", "xi2", "theta", "mu")
  } else if (frailty_type == "wlindley") {
    par_names <- c("zeta1", "delta1", "xi1", "zeta2", "delta2", "xi2", "s")
  } else {
    par_names <- c("zeta1", "delta1", "xi1", "zeta2", "delta2", "xi2", "theta")
  }
  if (n_cov > 0) par_names <- c(par_names, cov_names)
  n_pars <- length(par_names)
  n_pos_pars <- if (frailty_type == "genlindley") 8 else 7

  # Parameter transformation: unconstrained optimization
  unconstrained_to_natural <- function(u_pars) {
    nat_pars <- numeric(n_pars)
    if (frailty_type == "genlindley") {
      nat_pars[1:6] <- exp(u_pars[1:6])
      theta_val <- exp(u_pars[7])
      # mu is in (0, 1 + theta)
      mu_val <- (1 + theta_val) * stats::plogis(u_pars[8])
      nat_pars[7] <- theta_val
      nat_pars[8] <- mu_val
    } else {
      nat_pars[1:n_pos_pars] <- exp(u_pars[1:n_pos_pars])
    }
    if (n_cov > 0) {
      nat_pars[(n_pos_pars + 1):n_pars] <- u_pars[(n_pos_pars + 1):n_pars]
    }
    nat_pars
  }

  natural_to_unconstrained <- function(nat_pars) {
    u_pars <- numeric(n_pars)
    if (frailty_type == "genlindley") {
      u_pars[1:6] <- log(pmax(nat_pars[1:6], 1e-6))
      theta_val <- nat_pars[7]
      u_pars[7] <- log(pmax(theta_val, 1e-6))
      p_mu <- pmin(pmax(nat_pars[8] / (1 + theta_val), 1e-5), 1 - 1e-5)
      u_pars[8] <- stats::qlogis(p_mu)
    } else {
      u_pars[1:n_pos_pars] <- log(pmax(nat_pars[1:n_pos_pars], 1e-6))
    }
    if (n_cov > 0) {
      u_pars[(n_pos_pars + 1):n_pars] <- nat_pars[(n_pos_pars + 1):n_pars]
    }
    u_pars
  }

  # Default initial values
  if (is.null(init)) {
    if (baseline_dist == "gw") {
      base_init1 <- c(2.0, 0.05, 0.65)
      base_init2 <- c(2.7, 0.06, 0.67)
    } else {
      base_init1 <- c(3.8, 0.05, 1.2)
      base_init2 <- c(4.5, 0.06, 0.9)
    }
    if (frailty_type == "genlindley") frail_init <- c(2.3, 1.4)
    else if (frailty_type == "wlindley") frail_init <- c(1.1)
    else frail_init <- c(1.0)
    beta_init <- if (n_cov > 0) rep(0.01, n_cov) else numeric(0)
    init <- c(base_init1, base_init2, frail_init, beta_init)
  }

  init_u <- natural_to_unconstrained(init)

  obj_fn <- function(u_pars) {
    nat_pars <- unconstrained_to_natural(u_pars)
    ll <- loglik_add_frailty(nat_pars, t1, delta1, t2, delta2, X, frailty_type, baseline_dist)
    if (is.nan(ll) || is.na(ll) || ll <= -1e10) return(1e12)
    -ll
  }

  opt_res <- stats::optim(
    par = init_u,
    fn = obj_fn,
    method = method,
    hessian = hessian,
    control = control
  )

  est_nat <- unconstrained_to_natural(opt_res$par)
  names(est_nat) <- par_names
  max_ll <- -opt_res$value

  # Standard errors via Delta method / Jacobian
  vcov_mat <- matrix(NA_real_, nrow = n_pars, ncol = n_pars)
  rownames(vcov_mat) <- colnames(vcov_mat) <- par_names
  se_vec <- rep(NA_real_, n_pars)

  if (hessian && !is.null(opt_res$hessian)) {
    inv_hess <- tryCatch({
      solve(opt_res$hessian)
    }, error = function(e) {
      s_decomp <- svd(opt_res$hessian)
      d_vals <- s_decomp$d
      d_inv <- ifelse(d_vals > 1e-10 * max(d_vals), 1 / d_vals, 0)
      s_decomp$v %*% (d_inv * t(s_decomp$u))
    })

    # Jacobian matrix of transformation
    jac <- matrix(0, nrow = n_pars, ncol = n_pars)
    if (frailty_type == "genlindley") {
      for (k in 1:6) jac[k, k] <- est_nat[k]
      jac[7, 7] <- est_nat[7]
      p_mu <- stats::plogis(opt_res$par[8])
      jac[8, 7] <- p_mu
      jac[8, 8] <- (1 + est_nat[7]) * p_mu * (1 - p_mu)
    } else {
      for (k in 1:n_pos_pars) jac[k, k] <- est_nat[k]
    }
    if (n_cov > 0) {
      for (k in (n_pos_pars + 1):n_pars) jac[k, k] <- 1
    }

    vcov_mat <- jac %*% inv_hess %*% t(jac)
    rownames(vcov_mat) <- colnames(vcov_mat) <- par_names
    se_vec <- sqrt(pmax(0, diag(vcov_mat)))
  }

  z_vals <- est_nat / se_vec
  p_vals <- 2 * (1 - stats::pnorm(abs(z_vals)))
  ci_low <- est_nat - 1.96 * se_vec
  ci_high <- est_nat + 1.96 * se_vec

  # For strictly positive parameters, bound CI low at 0
  ci_low[1:n_pos_pars] <- pmax(0, ci_low[1:n_pos_pars])

  est_table <- data.frame(
    Estimate = est_nat,
    StdError = se_vec,
    z_value = z_vals,
    p_value = p_vals,
    CI_Lower = ci_low,
    CI_Upper = ci_high,
    row.names = par_names,
    stringsAsFactors = FALSE
  )

  # Criteria
  crit_list <- calc_model_criteria(
    loglik = max_ll,
    k = n_pars,
    n = length(t1)
  )

  # Kendall's tau and Frailty Variance
  if (frailty_type == "genlindley") {
    kt <- kendall_tau("genlindley", theta = est_nat[7], mu = est_nat[8])
    fv <- frailty_variance("genlindley", theta = est_nat[7], mu = est_nat[8])
  } else if (frailty_type == "wlindley") {
    kt <- kendall_tau("wlindley", s_param = est_nat[7])
    fv <- frailty_variance("wlindley", s_param = est_nat[7])
  } else if (frailty_type == "gamma") {
    kt <- kendall_tau("gamma", theta = est_nat[7])
    fv <- frailty_variance("gamma", theta = est_nat[7])
  } else {
    kt <- kendall_tau("invgaussian", theta = est_nat[7])
    fv <- frailty_variance("invgaussian", theta = est_nat[7])
  }

  out <- list(
    coefficients = est_nat,
    estimates = est_table,
    vcov = vcov_mat,
    loglik = max_ll,
    criteria = crit_list,
    kendall_tau = kt,
    frailty_variance = fv,
    convergence = opt_res$convergence,
    opt_output = opt_res,
    model_info = list(
      frailty_type = frailty_type,
      baseline_dist = baseline_dist,
      n_obs = length(t1),
      n_pars = n_pars,
      par_names = par_names,
      formula = formula,
      call = call
    )
  )

  class(out) <- c("add_frailty_mle", "add_frailty")
  out
}
