#' Fit Shared Gamma Frailty Model for k-Out-of-M One-Shot Devices
#'
#' @description
#' Fits the shared gamma frailty model with Exponential or Weibull component lifetimes
#' to one-shot device testing data with multiple failure modes under constant-stress
#' accelerated life tests. Computes maximum likelihood parameter estimates, standard
#' errors via the missing information principle, asymptotic confidence intervals,
#' and information criteria.
#'
#' @param formula Optional formula of the form \code{cbind(n_none, n_T, ...) ~ stress} or
#'   \code{cbind(n_none, n_T, ...) ~ stress + time}.
#' @param data A data frame containing the test group observations, covariates, inspection times,
#'   and malfunction subset counts.
#' @param time A character string specifying the column name in \code{data} containing inspection times \eqn{\tau_i}.
#' @param stress An optional character vector of column names in \code{data} containing stress covariates.
#' @param counts An optional character vector of \eqn{2^M} column names in \code{data} containing device counts
#'   for each subset in the power set \eqn{P(\Omega)}.
#' @param components An optional character vector of \eqn{M} component / failure mode names
#'   (e.g., \code{c("Turn", "Ground")}). If omitted, defaults to \code{c("Comp1", "Comp2", ...)}.
#' @param baseline A character string specifying the baseline lifetime distribution:
#'   \code{"exponential"} (default) or \code{"weibull"}.
#' @param frailty A logical value: \code{TRUE} (default) to fit the shared gamma frailty model (\eqn{\beta > 0}),
#'   or \code{FALSE} to fit the independence model (\eqn{\beta = 0}).
#' @param method A character string specifying the estimation algorithm: \code{"em"} (default, Expectation-Maximization)
#'   or \code{"mle"} (direct Quasi-Newton optimization).
#' @param start An optional numeric vector of starting values for the parameters.
#' @param ci_level A numeric scalar between 0 and 1 specifying the asymptotic confidence level (default 0.95).
#' @param control An optional list of control parameters passed to the optimizer (e.g. \code{maxit}, \code{tol}).
#' @param ... Additional arguments.
#'
#' @return An S3 object of class \code{"oneshot_kofm"} containing:
#'   \describe{
#'     \item{coefficients}{A named numeric vector of maximum likelihood estimates.}
#'     \item{vcov}{The asymptotic variance-covariance matrix of the estimators.}
#'     \item{se}{A named numeric vector of standard errors.}
#'     \item{confint}{A matrix with lower (\code{LCI}) and upper (\code{UCI}) asymptotic confidence bounds.}
#'     \item{loglik}{The maximized log-likelihood value.}
#'     \item{n_params}{The total number of estimated parameters.}
#'     \item{aic}{Akaike Information Criterion (\eqn{-2\log L + 2k}).}
#'     \item{bic}{Bayesian Information Criterion (\eqn{-2\log L + k \ln N}).}
#'     \item{aicc}{Corrected AIC for small sample sizes.}
#'     \item{hqic}{Hannan-Quinn Information Criterion.}
#'     \item{sample_size}{Total number of devices tested across all groups.}
#'     \item{fitted_probs}{An \eqn{I \times 2^M} matrix of model-estimated cell probabilities.}
#'     \item{fitted_counts}{An \eqn{I \times 2^M} matrix of expected device counts.}
#'     \item{residuals}{An \eqn{I \times 2^M} matrix of standardized Pearson residuals.}
#'     \item{kendall_tau}{A list with the estimated Kendall's tau association (\eqn{\tau = \beta / (\beta + 2)}),
#'       standard error, and confidence interval (if \code{frailty = TRUE}).}
#'     \item{M}{The number of components.}
#'     \item{components}{Character vector of component names.}
#'     \item{baseline}{Baseline distribution used (\code{"exponential"} or \code{"weibull"}).}
#'     \item{frailty}{Logical flag indicating whether gamma frailty was modeled.}
#'     \item{method}{Estimation method used (\code{"em"} or \code{"mle"}).}
#'     \item{data}{The processed input data frame.}
#'     \item{call}{The matched function call.}
#'   }
#'
#' @details
#' The shared gamma frailty model links the component hazard rates
#' \eqn{\lambda_{i,m} = \exp(a_m^T x_i)} with stress covariates.
#'
#' Missing values are checked upon input; observations with missing values in required
#' variables are rejected with an informative error message.
#'
#' @references
#' Ling, M.-H., Balakrishnan, N., Yu, C., and So, H. Y. (2021). Inference for One-Shot
#' Devices with Dependent k-Out-of-M Structured Components under Gamma Frailty.
#' \emph{Mathematics}, 9(23), 3032. \doi{10.3390/math9233032}.
#'
#' @examples
#' data(classH)
#' # Fit shared gamma frailty model with Exponential baseline
#' fit_exp <- oneshot_kofm(
#'   data = classH,
#'   time = "time",
#'   stress = "stress",
#'   counts = c("n_none", "n_T", "n_G", "n_both"),
#'   components = c("Turn", "Ground"),
#'   baseline = "exponential",
#'   frailty = TRUE,
#'   method = "em"
#' )
#' print(fit_exp)
#' summary(fit_exp)
#' @export
oneshot_kofm <- function(formula = NULL,
                         data,
                         time = "time",
                         stress = "stress",
                         counts = NULL,
                         components = NULL,
                         baseline = c("exponential", "weibull"),
                         frailty = TRUE,
                         method = c("em", "mle"),
                         start = NULL,
                         ci_level = 0.95,
                         control = list(),
                         ...) {
  cl <- match.call()
  baseline <- match.arg(baseline)
  method <- match.arg(method)
  
  if (missing(data) || is.null(data)) {
    stop("Argument 'data' must be supplied as a data frame.")
  }
  if (!is.data.frame(data)) {
    data <- as.data.frame(data)
  }
  
  if (is.null(frailty) || length(frailty) != 1 || is.na(frailty) || !is.logical(frailty)) {
    stop("'frailty' must be a single logical value (TRUE or FALSE).")
  }
  
  if (is.null(ci_level) || length(ci_level) != 1 || is.na(ci_level) || ci_level <= 0 || ci_level >= 1) {
    stop("'ci_level' must be a numeric scalar between 0 and 1 (e.g. 0.95).")
  }
  
  # Parse formula interface if provided
  if (!is.null(formula)) {
    mf <- stats::model.frame(formula, data = data)
    lhs <- stats::model.response(mf)
    
    if (is.null(lhs) || !is.matrix(lhs)) {
      resp_term <- formula[[2]]
      if (is.call(resp_term) && identical(as.character(resp_term[[1]]), "cbind")) {
        count_names <- as.character(as.list(resp_term)[-1])
        count_mat <- as.matrix(data[, count_names, drop = FALSE])
      } else {
        stop("The left-hand side of 'formula' must be a matrix of subset counts, e.g. cbind(n_none, n_T, n_G, n_both).")
      }
    } else {
      count_mat <- as.matrix(lhs)
    }
    
    mt <- stats::terms(formula, data = data)
    X_cov <- stats::model.matrix(mt, data = data)
    
    if (time %in% colnames(data)) {
      tau_vec <- as.numeric(data[[time]])
    } else if ("(time)" %in% colnames(mf)) {
      tau_vec <- as.numeric(mf[["(time)"]])
    } else {
      time_idx <- grep("time", colnames(X_cov), ignore.case = TRUE)
      if (length(time_idx) > 0) {
        tau_vec <- as.numeric(X_cov[, time_idx[1]])
        X_cov <- X_cov[, -time_idx[1], drop = FALSE]
      } else {
        stop("Inspection time column '", time, "' was not found in 'data'.")
      }
    }
  } else {
    # Argument interface
    if (!(time %in% colnames(data))) {
      stop("Inspection time column '", time, "' was not found in 'data'.")
    }
    tau_vec <- as.numeric(data[[time]])
    
    if (!is.null(stress)) {
      missing_stress <- setdiff(stress, colnames(data))
      if (length(missing_stress) > 0) {
        stop("Stress covariate column(s) not found in 'data': ", paste(missing_stress, collapse = ", "))
      }
      X_cov <- cbind(1, as.matrix(data[, stress, drop = FALSE]))
      colnames(X_cov) <- c("(Intercept)", stress)
    } else {
      X_cov <- matrix(1, nrow = nrow(data), ncol = 1)
      colnames(X_cov) <- "(Intercept)"
    }
    
    if (is.null(counts)) {
      numeric_cols <- vapply(data, is.numeric, logical(1))
      cand_cols <- names(numeric_cols)[numeric_cols]
      cand_cols <- setdiff(cand_cols, c(time, stress))
      n_cand <- length(cand_cols)
      
      M_guess <- log2(n_cand)
      if (M_guess == as.integer(M_guess) && M_guess >= 1) {
        counts <- cand_cols
      } else {
        stop("Please specify 'counts' as a character vector of 2^M subset count column names.")
      }
    }
    
    missing_counts <- setdiff(counts, colnames(data))
    if (length(missing_counts) > 0) {
      stop("Count column(s) not found in 'data': ", paste(missing_counts, collapse = ", "))
    }
    count_mat <- as.matrix(data[, counts, drop = FALSE])
  }
  
  if (any(is.na(count_mat)) || any(count_mat < 0)) {
    stop("'count_mat' must contain non-negative numeric counts without missing values.")
  }
  if (any(is.na(X_cov))) {
    stop("Covariates 'X_cov' must not contain missing values.")
  }
  if (any(is.na(tau_vec)) || any(tau_vec < 0)) {
    stop("Inspection times 'tau' must be non-negative without missing values.")
  }
  
  n_subsets <- ncol(count_mat)
  M <- as.integer(round(log2(n_subsets)))
  if (2^M != n_subsets || M < 1) {
    stop("Number of count columns (", n_subsets, ") is not a power of 2 (2^M).")
  }
  
  if (is.null(components)) {
    components <- paste0("Comp", seq_len(M))
  } else if (length(components) != M) {
    stop("Length of 'components' (", length(components), ") must equal M (", M, ").")
  }
  
  p_cov <- ncol(X_cov)
  I_groups <- nrow(X_cov)
  powerset <- generate_powerset(M)
  
  cov_names <- colnames(X_cov)
  if (is.null(cov_names)) cov_names <- paste0("x", seq_len(p_cov) - 1)
  
  par_names <- character(0)
  for (m in seq_len(M)) {
    par_names <- c(par_names, paste0("a_", components[m], "_", cov_names))
  }
  if (baseline == "weibull") {
    for (m in seq_len(M)) {
      par_names <- c(par_names, paste0("b_", components[m], "_", cov_names))
    }
  }
  if (frailty) {
    par_names <- c(par_names, "beta")
  }
  
  max_iter <- if (!is.null(control$maxit)) control$maxit else 500
  tol <- if (!is.null(control$tol)) control$tol else 1e-5
  trace <- if (!is.null(control$trace)) control$trace else FALSE
  
  if (method == "em") {
    fit_res <- em_oneshot_kofm(
      count_mat = count_mat,
      X_cov = X_cov,
      tau = tau_vec,
      M = M,
      baseline = baseline,
      frailty = frailty,
      start = start,
      max_iter = max_iter,
      tol = tol,
      trace = trace
    )
    est_par <- fit_res$par
    max_loglik <- fit_res$loglik
  } else {
    if (is.null(start)) {
      a_init <- matrix(0, nrow = M, ncol = p_cov)
      for (m in seq_len(M)) {
        a_init[m, 1] <- -log(mean(tau_vec) + 1)
        if (p_cov > 1) a_init[m, 2:p_cov] <- 0.01
      }
      b_init <- if (baseline == "weibull") matrix(0, nrow = M, ncol = p_cov) else NULL
      beta_init <- if (frailty) 0.3 else numeric(0)
      start <- c(as.vector(t(a_init)),
                 if (baseline == "weibull") as.vector(t(b_init)) else numeric(0),
                 beta_init)
    }
    
    n_par <- length(start)
    obj_fun <- function(p) {
      -loglik_oneshot(p, count_mat, X_cov, tau_vec, M, baseline, frailty, powerset)
    }
    
    if (frailty) {
      lower_bounds <- rep(-Inf, n_par)
      lower_bounds[n_par] <- 1e-5
      upper_bounds <- rep(Inf, n_par)
      upper_bounds[n_par] <- 0.5
      
      opt_step <- tryCatch({
        stats::optim(start, obj_fun, method = "L-BFGS-B",
                     lower = lower_bounds, upper = upper_bounds,
                     control = list(maxit = max_iter, factr = 1e7, pgtol = 1e-8))
      }, error = function(e) {
        stats::optim(start, obj_fun, method = "Nelder-Mead",
                     control = list(maxit = max_iter))
      })
    } else {
      opt_step <- stats::optim(start, obj_fun, method = "BFGS",
                               control = list(maxit = max_iter, reltol = 1e-8))
    }
    est_par <- opt_step$par
    max_loglik <- -opt_step$value
  }
  
  names(est_par) <- par_names
  n_params <- length(est_par)
  
  # Information matrix and standard errors
  info_mat <- fisher_info_oneshot(est_par, count_mat, X_cov, tau_vec, M, baseline, frailty, powerset)
  dimnames(info_mat) <- list(par_names, par_names)
  
  vcov_mat <- tryCatch({
    solve(info_mat)
  }, error = function(e) {
    MASS_ginv <- function(m) {
      s <- svd(m)
      pos <- s$d > max(1e-8 * s$d[1], 0)
      if (all(pos)) s$v %*% (1 / s$d * t(s$u)) else s$v[, pos, drop = FALSE] %*% ((1 / s$d[pos]) * t(s$u[, pos, drop = FALSE]))
    }
    MASS_ginv(info_mat)
  })
  dimnames(vcov_mat) <- list(par_names, par_names)
  
  se_vec <- sqrt(pmax(diag(vcov_mat), 0))
  names(se_vec) <- par_names
  
  # Confidence intervals
  z_crit <- stats::qnorm(1 - (1 - ci_level) / 2)
  lci <- est_par - z_crit * se_vec
  uci <- est_par + z_crit * se_vec
  
  if (frailty) {
    lci["beta"] <- max(0, lci["beta"])
    uci["beta"] <- min(0.5, max(0, uci["beta"]))
  }
  confint_mat <- cbind(LCI = lci, UCI = uci)
  rownames(confint_mat) <- par_names
  
  # Sample size and information criteria
  K_i <- rowSums(count_mat)
  N_tot <- sum(K_i)
  
  aic_val <- -2 * max_loglik + 2 * n_params
  bic_val <- -2 * max_loglik + n_params * log(N_tot)
  aicc_val <- if (N_tot > n_params + 1) aic_val + (2 * n_params * (n_params + 1)) / (N_tot - n_params - 1) else aic_val
  hqic_val <- -2 * max_loglik + 2 * n_params * log(log(N_tot))
  
  # Fitted probabilities and counts
  prob_res <- compute_probabilities(est_par, X_cov, tau_vec, M, baseline, frailty, powerset)
  fitted_probs <- prob_res$prob_matrix
  fitted_counts <- fitted_probs * K_i
  
  # Pearson residuals
  residuals_mat <- (count_mat - fitted_counts) / sqrt(pmax(fitted_counts, 1e-8))
  
  # Kendall's tau association list
  kendall_tau <- if (frailty && "beta" %in% names(est_par)) {
    beta_est <- est_par["beta"]
    tau_est <- beta_est / (beta_est + 2)
    # Delta method SE for tau
    # d(tau)/d(beta) = 2 / (beta + 2)^2
    d_tau <- 2 / (beta_est + 2)^2
    se_tau <- as.numeric(d_tau * se_vec["beta"])
    list(
      estimate = tau_est,
      se = se_tau,
      confint = c(LCI = max(0, tau_est - z_crit * se_tau),
                  UCI = min(0.2, tau_est + z_crit * se_tau))
    )
  } else {
    NULL
  }
  
  structure(
    list(
      coefficients = est_par,
      vcov = vcov_mat,
      se = se_vec,
      confint = confint_mat,
      loglik = max_loglik,
      n_params = n_params,
      aic = aic_val,
      bic = bic_val,
      aicc = aicc_val,
      hqic = hqic_val,
      sample_size = N_tot,
      fitted_probs = fitted_probs,
      fitted_counts = fitted_counts,
      residuals = residuals_mat,
      kendall_tau = kendall_tau,
      M = M,
      components = components,
      baseline = baseline,
      frailty = frailty,
      method = method,
      ci_level = ci_level,
      data = data,
      call = cl
    ),
    class = "oneshot_kofm"
  )
}
