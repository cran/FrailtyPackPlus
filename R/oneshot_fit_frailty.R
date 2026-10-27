#' Fit Shared Gamma Frailty Model for One-Shot Device Testing Data
#'
#' @description
#' Fits the Weibull or Exponential shared gamma frailty model to one-shot device
#' testing data with multiple failure modes under constant-stress accelerated life tests,
#' computing maximum likelihood parameter estimates, standard errors, asymptotic
#' confidence intervals, and model selection criteria.
#'
#' @param formula Optional formula of the form \code{cbind(n_none, n_T, ...) ~ stress} or
#'   \code{cbind(n_none, n_T, ...) ~ stress + time}.
#' @param data A data frame containing the test group observations, covariates, inspection times,
#'   and malfunction subset counts.
#' @param time A character string specifying the column name in \code{data} containing inspection times \eqn{\tau_i}.
#' @param stress An optional character vector of column names in \code{data} containing stress covariates.
#' @param counts An optional character vector of \eqn{2^M} column names in \code{data} containing device counts
#'   for each subset in the power set \eqn{P(Omega)}.
#' @param components An optional character vector of \eqn{M} component / failure mode names
#'   (e.g., \code{c("Turn", "Ground")}). If omitted, defaults to \code{c("Comp1", "Comp2", ...)}.
#' @param baseline A character string specifying the baseline lifetime distribution:
#'   \code{"weibull"} (default) or \code{"exponential"}.
#' @param frailty A logical value: \code{TRUE} (default) to fit the shared gamma frailty model (\eqn{\beta > 0}),
#'   or \code{FALSE} to fit the independence model (\eqn{\beta = 0}).
#' @param method A character string specifying the estimation algorithm: \code{"mle"} (default, direct Quasi-Newton)
#'   or \code{"em"} (Expectation-Maximization algorithm).
#' @param start An optional numeric vector of starting values for the parameters.
#' @param ci_level A numeric scalar between 0 and 1 specifying the asymptotic confidence level (default 0.95).
#' @param control An optional list of control parameters passed to the optimizer (e.g. \code{maxit}, \code{tol}).
#' @param ... Additional arguments.
#'
#' @return An S3 object of class \code{"oneshot_frailty"} containing:
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
#'     \item{fitted_probs}{An \eqn{I x 2^M} matrix of model-estimated cell probabilities.}
#'     \item{fitted_counts}{An \eqn{I x 2^M} matrix of expected device counts.}
#'     \item{residuals}{An \eqn{I x 2^M} matrix of standardized Pearson residuals.}
#'     \item{kendall_tau}{A list with the estimated Kendall's tau association (\eqn{\tau = \beta / (\beta + 2)}),
#'       standard error, and confidence interval (if \code{frailty = TRUE}).}
#'     \item{M}{The number of components.}
#'     \item{components}{Character vector of component names.}
#'     \item{baseline}{Baseline distribution used (\code{"weibull"} or \code{"exponential"}).}
#'     \item{frailty}{Logical flag indicating whether gamma frailty was modeled.}
#'     \item{method}{Estimation method used (\code{"mle"} or \code{"em"}).}
#'     \item{data}{The processed input data frame.}
#'     \item{call}{The matched function call.}
#'   }
#'
#' @details
#' The shared gamma frailty model links the component baseline hazard functions
#' \eqn{h_{0,i,m}(t) = \alpha_{i,m} \lambda_{i,m} (\lambda_{i,m} t)^{\alpha_{i,m}-1}}
#' with stress covariates via log-linear links \eqn{\lambda_{i,m} = \exp(a_m^T x_i)}
#' and \eqn{\alpha_{i,m} = \exp(b_m^T x_i)}.
#'
#' Missing values are checked upon input; observations with missing values in required
#' variables are rejected with an informative error message.
#'
#' @references
#' Yu, C., So, H. Y., and Balakrishnan, N. (2025). \emph{Inference for one-shot devices
#' with Weibull component lifetimes under gamma frailty model}. In: Stochastic Modeling
#' and Statistical Methods, Elsevier Inc., Chapter 1, pp. 1--34.
#' \doi{10.1016/B978-0-44-331694-4.00006-2}.
#'
#' @examples
#' data(classH)
#' # Fit shared gamma frailty model with Weibull baseline
#' fit_w <- oneshot_frailty(
#'   data = classH,
#'   time = "time",
#'   stress = "stress",
#'   counts = c("n_none", "n_T", "n_G", "n_both"),
#'   components = c("Turn", "Ground"),
#'   baseline = "weibull",
#'   frailty = TRUE
#' )
#' print(fit_w)
#' summary(fit_w)
#' @export
oneshot_frailty <- function(formula = NULL,
                            data,
                            time = "time",
                            stress = "stress",
                            counts = NULL,
                            components = NULL,
                            baseline = c("weibull", "exponential"),
                            frailty = TRUE,
                            method = c("mle", "em"),
                            start = NULL,
                            ci_level = 0.95,
                            control = list(),
                            ...) {
  cl <- match.call()
  baseline <- match.arg(baseline)
  method <- match.arg(method)
  
  if (missing(data) || is.null(data) || !is.data.frame(data)) {
    stop("'data' must be a valid non-empty data frame.")
  }
  
  if (!is.null(formula)) {
    mf <- stats::model.frame(formula, data = data)
    y_mat <- stats::model.response(mf)
    if (!is.matrix(y_mat)) {
      stop("LHS of formula must be a matrix of counts created with cbind().")
    }
    count_mat <- y_mat
    cov_terms <- stats::terms(formula)
    X_cov <- stats::model.matrix(cov_terms, data = mf)
    
    if (!(time %in% colnames(data))) {
      stop(sprintf("Time variable '%s' not found in data.", time))
    }
    tau <- as.numeric(data[[time]])
  } else {
    if (!(time %in% colnames(data))) {
      stop(sprintf("Time variable '%s' not found in data.", time))
    }
    tau <- as.numeric(data[[time]])
    
    if (is.null(stress) || length(stress) == 0) {
      X_cov <- matrix(1, nrow = nrow(data), ncol = 1)
      colnames(X_cov) <- "(Intercept)"
    } else {
      missing_vars <- setdiff(stress, colnames(data))
      if (length(missing_vars) > 0) {
        stop(sprintf("Covariate(s) not found in data: %s", paste(missing_vars, collapse = ", ")))
      }
      cov_df <- data[, stress, drop = FALSE]
      X_cov <- stats::model.matrix(~ ., data = cov_df)
    }
    
    if (is.null(counts)) {
      count_cols <- grep("^n_", colnames(data), value = TRUE)
      if (length(count_cols) == 0) {
        stop("Please specify 'counts' column names representing the 2^M subset counts.")
      }
      count_mat <- as.matrix(data[, count_cols, drop = FALSE])
    } else {
      missing_counts <- setdiff(counts, colnames(data))
      if (length(missing_counts) > 0) {
        stop(sprintf("Count column(s) not found in data: %s", paste(missing_counts, collapse = ", ")))
      }
      count_mat <- as.matrix(data[, counts, drop = FALSE])
    }
  }
  
  if (any(is.na(X_cov)) || any(is.na(tau)) || any(is.na(count_mat))) {
    stop("Missing values (NA/NaN) detected in data or covariates. Please clean data before fitting.")
  }
  if (any(tau <= 0)) {
    stop("Inspection times 'tau' must be strictly positive.")
  }
  if (any(count_mat < 0)) {
    stop("Observed counts must be non-negative integers.")
  }
  
  n_subsets <- ncol(count_mat)
  M <- round(log2(n_subsets))
  if (2^M != n_subsets || M < 1) {
    stop(sprintf("The number of count columns (%d) must be equal to 2^M for some integer M >= 1.", n_subsets))
  }
  
  if (is.null(components)) {
    components <- paste0("Comp", seq_len(M))
  } else if (length(components) != M) {
    stop(sprintf("'components' must have length %d (matching M).", M))
  }
  
  powerset <- generate_powerset(M)
  p_cov <- ncol(X_cov)
  cov_names <- colnames(X_cov)
  if (is.null(cov_names)) {
    cov_names <- paste0("x", seq_len(p_cov) - 1)
  }
  
  param_names <- character(0)
  for (m in seq_len(M)) {
    for (p in seq_len(p_cov)) {
      param_names <- c(param_names, sprintf("a_%s_%s", components[m], cov_names[p]))
    }
  }
  if (baseline == "weibull") {
    for (m in seq_len(M)) {
      for (p in seq_len(p_cov)) {
        param_names <- c(param_names, sprintf("b_%s_%s", components[m], cov_names[p]))
      }
    }
  }
  if (frailty) {
    param_names <- c(param_names, "beta")
  }
  
  n_params <- length(param_names)
  
  obj_fun <- function(p) {
    -loglik_oneshot(p, count_mat, X_cov, tau, M, baseline, frailty, powerset)
  }
  
  # Smart multi-start heuristic if start is NULL
  if (is.null(start)) {
    mean_log_tau <- log(mean(tau))
    candidate_starts <- list()
    
    # Candidate 1: centered around standard rates
    c1 <- rep(0, n_params)
    for (m in seq_len(M)) {
      c1[(m - 1) * p_cov + 1] <- -mean_log_tau
    }
    if (frailty) c1[n_params] <- 0.1
    candidate_starts[[1]] <- c1
    
    # Candidate 2: small negative slope if stress present
    if (p_cov > 1) {
      c2 <- c1
      for (m in seq_len(M)) {
        c2[(m - 1) * p_cov + 1] <- -mean_log_tau * 1.2
        c2[(m - 1) * p_cov + 2] <- 0.01
      }
      if (frailty) c2[n_params] <- 0.08
      candidate_starts[[2]] <- c2
      
      c3 <- c1
      for (m in seq_len(M)) {
        c3[(m - 1) * p_cov + 1] <- -mean_log_tau * 0.8
        c3[(m - 1) * p_cov + 2] <- -0.01
      }
      if (frailty) c3[n_params] <- 0.2
      candidate_starts[[3]] <- c3
    }
    
    best_init <- candidate_starts[[1]]
    best_v <- obj_fun(best_init)
    for (cand in candidate_starts) {
      v_curr <- obj_fun(cand)
      if (v_curr < best_v) {
        best_v <- v_curr
        best_init <- cand
      }
    }
    start <- best_init
  }
  names(start) <- param_names
  
  if (method == "em") {
    max_it <- if (!is.null(control$maxit)) control$maxit else 500
    tol_val <- if (!is.null(control$tol)) control$tol else 1e-5
    em_res <- em_oneshot_frailty(count_mat, X_cov, tau, M, baseline, frailty,
                                 start = start, max_iter = max_it, tol = tol_val)
    est_par <- em_res$par
    max_ll <- em_res$loglik
  } else {
    if (frailty) {
      lower_bounds <- rep(-Inf, n_params)
      lower_bounds[n_params] <- 1e-6
      upper_bounds <- rep(Inf, n_params)
      
      opt_res <- tryCatch({
        stats::optim(start, obj_fun, method = "L-BFGS-B",
                     lower = lower_bounds, upper = upper_bounds,
                     control = list(maxit = 1000, factr = 1e7, pgtol = 1e-8))
      }, error = function(e) {
        stats::optim(start, obj_fun, method = "Nelder-Mead",
                     control = list(maxit = 2000))
      })
    } else {
      opt_res <- stats::optim(start, obj_fun, method = "BFGS",
                              control = list(maxit = 1000, reltol = 1e-10))
    }
    est_par <- opt_res$par
    max_ll <- -opt_res$value
  }
  
  names(est_par) <- param_names
  
  nll_wrapper <- function(p) {
    -loglik_oneshot(p, count_mat, X_cov, tau, M, baseline, frailty, powerset)
  }
  
  hess <- tryCatch({
    numDeriv::hessian(nll_wrapper, est_par)
  }, error = function(e) {
    matrix(NA, n_params, n_params)
  })
  
  vcov_mat <- tryCatch({
    inv_h <- solve(hess)
    inv_h <- (inv_h + t(inv_h)) / 2
    inv_h
  }, error = function(e) {
    diag(1e-4, n_params)
  })
  
  dimnames(vcov_mat) <- list(param_names, param_names)
  
  var_diag <- diag(vcov_mat)
  var_diag[is.na(var_diag) | var_diag < 0] <- 1e-4
  se_vec <- sqrt(var_diag)
  names(se_vec) <- param_names
  
  z_crit <- stats::qnorm(1 - (1 - ci_level) / 2)
  ci_mat <- cbind(
    LCI = est_par - z_crit * se_vec,
    UCI = est_par + z_crit * se_vec
  )
  if (frailty) {
    ci_mat["beta", "LCI"] <- max(0, ci_mat["beta", "LCI"])
  }
  rownames(ci_mat) <- param_names
  
  total_devices <- sum(count_mat)
  aic_val <- -2 * max_ll + 2 * n_params
  bic_val <- -2 * max_ll + n_params * log(total_devices)
  aicc_val <- if (total_devices > n_params + 1) {
    aic_val + (2 * n_params * (n_params + 1)) / (total_devices - n_params - 1)
  } else {
    aic_val
  }
  hqic_val <- -2 * max_ll + 2 * n_params * log(log(max(total_devices, 3)))
  
  fit_probs <- compute_probabilities(est_par, X_cov, tau, M, baseline, frailty, powerset)$prob_matrix
  row_totals <- rowSums(count_mat)
  fit_counts <- fit_probs * matrix(rep(row_totals, n_subsets), nrow = nrow(count_mat), ncol = n_subsets)
  
  denom <- sqrt(pmax(fit_counts * (1 - fit_probs), 1e-10))
  pearson_res <- (count_mat - fit_counts) / denom
  
  subset_labels <- vapply(powerset, function(s) {
    if (length(s) == 0) "none"
    else paste(components[s], collapse = "+")
  }, character(1))
  colnames(fit_probs) <- subset_labels
  colnames(fit_counts) <- subset_labels
  colnames(pearson_res) <- subset_labels
  
  if (frailty) {
    beta_est <- est_par["beta"]
    tau_est <- beta_est / (beta_est + 2)
    se_beta <- se_vec["beta"]
    se_tau <- se_beta * (2 / (beta_est + 2)^2)
    tau_ci <- c(
      LCI = max(0, tau_est - z_crit * se_tau),
      UCI = min(1, tau_est + z_crit * se_tau)
    )
    kendall_info <- list(tau = tau_est, se = se_tau, confint = tau_ci)
  } else {
    kendall_info <- list(tau = 0, se = 0, confint = c(LCI = 0, UCI = 0))
  }
  
  out <- list(
    coefficients = est_par,
    vcov = vcov_mat,
    se = se_vec,
    confint = ci_mat,
    loglik = max_ll,
    n_params = n_params,
    aic = aic_val,
    bic = bic_val,
    aicc = aicc_val,
    hqic = hqic_val,
    sample_size = total_devices,
    fitted_probs = fit_probs,
    fitted_counts = fit_counts,
    residuals = pearson_res,
    kendall_tau = kendall_info,
    M = M,
    components = components,
    baseline = baseline,
    frailty = frailty,
    method = method,
    data = data,
    X_cov = X_cov,
    tau = tau,
    count_mat = count_mat,
    ci_level = ci_level,
    call = cl
  )
  class(out) <- "oneshot_frailty"
  out
}
