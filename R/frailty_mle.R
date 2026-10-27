#' Maximum Likelihood Estimation of Shared Frailty Models
#'
#' Estimates shared, additive, multiplicative, and merged frailty models
#' using maximum likelihood estimation with analytical/numerical Hessian.
#'
#' @param formula Formula expression for bivariate or univariate survival data.
#' @param data Data frame containing survival data and covariates.
#' @param baseline Baseline distribution name.
#' @param frailty Frailty distribution name.
#' @param type Model formulation type.
#' @param init Initial parameter vector (optional).
#' @param method Optimization method passed to \code{stats::optim} (default \code{"BFGS"}).
#' @param hessian Logical; if \code{TRUE}, computes the Hessian matrix.
#' @param ... Additional arguments.
#'
#' @return An S3 object of class \code{"frailty_mle"}.
#' @export
#' @examples
#' data(kidney_frailty)
#' fit_mle <- frailty_mle(
#'   Surv(time1, time2, status1, status2) ~ age + sex,
#'   data = kidney_frailty,
#'   baseline = "weibull",
#'   frailty = "gamma",
#'   type = "multiplicative_bivariate"
#' )
#' print(fit_mle)
#' summary(fit_mle)
frailty_mle <- function(formula, data, baseline = c("weibull", "exponential",
                                                   "genexponential", "exppower",
                                                   "lfr", "genloglogistic",
                                                   "genweibull", "genrayleigh",
                                                   "weightedexp", "extendedweibull"),
                        frailty = c("gamma", "invgauss", "posstable",
                                    "compoundpoisson", "compoundnb",
                                    "genlindley", "weightedlindley",
                                    "mixture_invgauss"),
                        type = c("multiplicative_bivariate", "additive_bivariate",
                                 "multiplicative_univariate", "additive_univariate",
                                 "merged_univariate"),
                        init = NULL, method = "BFGS", hessian = TRUE, ...) {
  cl <- match.call()
  baseline <- match.arg(baseline)
  frailty <- match.arg(frailty)
  type <- match.arg(type)
  
  parsed <- parse_frailty_data(formula, data, type)
  n_base <- get_n_base_params(baseline)
  n_frail <- get_n_frail_params(frailty)
  n_cov <- parsed$p_cov
  n_total <- n_base + n_cov + n_frail
  
  param_names <- get_param_names(baseline, frailty, parsed$covar_names)
  
  if (is.null(init)) {
    init_base <- rep(1.0, n_base)
    if (baseline == "weibull") init_base <- c(0.01, 1.0)
    if (baseline == "exppower") init_base <- c(0.01, 1.0)
    init_beta <- rep(0.0, n_cov)
    init_frail <- rep(1.0, n_frail)
    if (frailty == "posstable") init_frail <- 0.7
    if (frailty == "mixture_invgauss") init_frail <- c(0.5, 1.0, 1.0)
    init <- c(init_base, init_beta, init_frail)
  }
  
  trans_info <- get_param_transforms(baseline, frailty, n_cov)
  to_uncon <- trans_info$to_uncon
  to_con <- trans_info$to_con
  init_uncon <- to_uncon(init)
  
  # Negative log-likelihood objective on unconstrained space
  neg_loglik_fn <- function(uncon_p) {
    con_p <- to_con(uncon_p)
    ll <- loglik_frailty(con_p, parsed, baseline = baseline,
                         frailty = frailty, type = type, return_indiv = FALSE)
    if (!is.finite(ll)) return(1e10)
    -ll
  }
  
  opt <- stats::optim(init_uncon, neg_loglik_fn, method = method, hessian = hessian, ...)
  
  est_con <- to_con(opt$par)
  names(est_con) <- param_names
  
  # Variance-covariance matrix
  vcov_mat <- matrix(NA, nrow = n_total, ncol = n_total,
                     dimnames = list(param_names, param_names))
  se_vec <- rep(NA, n_total)
  
  if (hessian && !is.null(opt$hessian)) {
    h <- opt$hessian
    eig <- eigen(h, symmetric = TRUE)
    if (all(eig$values > 1e-8)) {
      inv_h_uncon <- solve(h)
      # Delta method for Jacobian of transformation
      # For positive param: con = exp(uncon), dcon/duncon = con
      # For logit: con = plogis(uncon), dcon/duncon = con*(1-con)
      # For beta: 1
      J_diag <- numeric(n_total)
      J_diag[1:n_base] <- est_con[1:n_base]
      if (n_cov > 0) J_diag[(n_base + 1):(n_base + n_cov)] <- 1
      idx_f <- (n_base + n_cov + 1):(n_base + n_cov + n_frail)
      if (frailty == "posstable") {
        a <- est_con[idx_f[1]]
        J_diag[idx_f[1]] <- a * (1 - a)
      } else if (frailty == "mixture_invgauss") {
        p1 <- est_con[idx_f[1]]
        J_diag[idx_f[1]] <- p1 * (1 - p1)
        J_diag[idx_f[2]] <- est_con[idx_f[2]]
        J_diag[idx_f[3]] <- est_con[idx_f[3]]
      } else {
        J_diag[idx_f] <- est_con[idx_f]
      }
      
      J_mat <- diag(J_diag, n_total)
      vcov_mat <- J_mat %*% inv_h_uncon %*% J_mat
      dimnames(vcov_mat) <- list(param_names, param_names)
      se_vec <- sqrt(pmax(diag(vcov_mat), 0))
    }
  }
  
  log_lik_val <- -opt$value
  aic_val <- 2 * n_total - 2 * log_lik_val
  bic_val <- n_total * log(parsed$n) - 2 * log_lik_val
  
  summary_df <- data.frame(
    Estimate = est_con,
    Std_Error = se_vec,
    z_value = ifelse(!is.na(se_vec) & se_vec > 0, est_con / se_vec, NA),
    p_value = ifelse(!is.na(se_vec) & se_vec > 0, 2 * (1 - stats::pnorm(abs(est_con / se_vec))), NA),
    CI_Lower = ifelse(!is.na(se_vec), est_con - 1.96 * se_vec, NA),
    CI_Upper = ifelse(!is.na(se_vec), est_con + 1.96 * se_vec, NA),
    stringsAsFactors = FALSE
  )
  rownames(summary_df) <- param_names
  
  out <- list(
    coefficients = est_con,
    vcov = vcov_mat,
    se = se_vec,
    summary = summary_df,
    logLik = log_lik_val,
    AIC = aic_val,
    BIC = bic_val,
    convergence = opt$convergence,
    counts = opt$counts,
    message = opt$message,
    data = parsed,
    call = cl,
    baseline = baseline,
    frailty = frailty,
    type = type,
    param_names = param_names,
    n_params = n_total,
    n_obs = parsed$n
  )
  class(out) <- "frailty_mle"
  out
}
