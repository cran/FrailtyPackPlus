#' Maximum Likelihood Estimation Engine for Correlated Frailty Models
#'
#' @description
#' Maximum likelihood estimation via numeric optimization for parametric
#' correlated frailty bivariate survival models.
#'
#' @name mle_engine
#' @keywords internal
NULL

fit_corrFrailty_mle <- function(t1, d1, t2, d2, X,
                                baseline1, baseline2, frailty,
                                control) {
  n_obs <- length(t1)
  p_cov <- if (is.null(X)) 0 else ncol(X)
  cov_names <- if (p_cov > 0) colnames(X) else character(0)
  
  bnames1 <- paste0(baseline_param_names(baseline1), "_1")
  bnames2 <- paste0(baseline_param_names(baseline2), "_2")
  fnames <- frailty_param_names(frailty)
  all_param_names <- c(bnames1, bnames2, fnames, cov_names)
  n_params <- length(all_param_names)
  
  unpack_params <- function(theta) {
    nb1 <- length(bnames1); nb2 <- length(bnames2); nf <- length(fnames)
    idx <- 1
    bp1 <- theta[idx:(idx + nb1 - 1)]; idx <- idx + nb1
    bp2 <- theta[idx:(idx + nb2 - 1)]; idx <- idx + nb2
    fp <- if (nf > 0) {
      res <- theta[idx:(idx + nf - 1)]; idx <- idx + nf; res
    } else {
      numeric(0)
    }
    beta <- if (p_cov > 0) theta[idx:(idx + p_cov - 1)] else numeric(0)
    list(base1 = bp1, base2 = bp2, frailty = fp, beta = beta)
  }
  
  calc_loglik_indiv <- function(unpacked) {
    if (p_cov > 0) {
      eta_lin <- as.vector(X %*% unpacked$beta)
      eta1 <- exp(pmin(pmax(eta_lin, -20), 20))
      eta2 <- eta1
    } else {
      eta1 <- rep(1, n_obs)
      eta2 <- rep(1, n_obs)
    }
    bivariate_loglik_contrib(t1, d1, t2, d2, eta1, eta2,
                             baseline1, unpacked$base1,
                             baseline2, unpacked$base2,
                             frailty, unpacked$frailty)
  }
  
  to_unconstrained <- function(theta) {
    unp <- unpack_params(theta)
    nb1 <- length(bnames1); nb2 <- length(bnames2)
    phi <- numeric(n_params)
    idx <- 1
    phi[idx:(idx + nb1 - 1)] <- log(pmax(1e-6, unp$base1)); idx <- idx + nb1
    phi[idx:(idx + nb2 - 1)] <- log(pmax(1e-6, unp$base2)); idx <- idx + nb2
    if (frailty %in% c("gamma", "invgaussian")) {
      phi[idx] <- log(pmax(1e-6, unp$frailty[1]))
      phi[idx + 1] <- log(unp$frailty[2] / (1 - unp$frailty[2]))
      idx <- idx + 2
    } else if (frailty == "posstable") {
      phi[idx] <- log(unp$frailty[1] / (1 - unp$frailty[1]))
      phi[idx + 1] <- log(unp$frailty[2] / (1 - unp$frailty[2]))
      idx <- idx + 2
    } else if (frailty == "comppoisson") {
      phi[idx:(idx + 2)] <- log(pmax(1e-6, unp$frailty))
      idx <- idx + 3
    }
    if (p_cov > 0) {
      phi[idx:(idx + p_cov - 1)] <- unp$beta
    }
    phi
  }
  
  to_constrained <- function(phi) {
    nb1 <- length(bnames1); nb2 <- length(bnames2)
    theta <- numeric(n_params)
    idx <- 1
    theta[idx:(idx + nb1 - 1)] <- exp(phi[idx:(idx + nb1 - 1)]); idx <- idx + nb1
    theta[idx:(idx + nb2 - 1)] <- exp(phi[idx:(idx + nb2 - 1)]); idx <- idx + nb2
    if (frailty %in% c("gamma", "invgaussian")) {
      theta[idx] <- exp(phi[idx])
      theta[idx + 1] <- 1 / (1 + exp(-phi[idx + 1]))
      idx <- idx + 2
    } else if (frailty == "posstable") {
      theta[idx] <- 1 / (1 + exp(-phi[idx]))
      theta[idx + 1] <- 1 / (1 + exp(-phi[idx + 1]))
      idx <- idx + 2
    } else if (frailty == "comppoisson") {
      theta[idx:(idx + 2)] <- exp(phi[idx:(idx + 2)])
      idx <- idx + 3
    }
    if (p_cov > 0) {
      theta[idx:(idx + p_cov - 1)] <- phi[idx:(idx + p_cov - 1)]
    }
    theta
  }
  
  neg_loglik_obj <- function(phi) {
    theta <- to_constrained(phi)
    unp <- unpack_params(theta)
    ll_vec <- calc_loglik_indiv(unp)
    if (any(is.nan(ll_vec)) || any(is.na(ll_vec))) return(1e10)
    -sum(ll_vec)
  }
  
  # Default initial values
  rate1 <- max(0.01, sum(d1) / sum(t1))
  rate2 <- max(0.01, sum(d2) / sum(t2))
  init_bp1 <- rep(1.0, length(bnames1))
  init_bp1[which(grepl("lambda", bnames1))] <- rate1
  init_bp2 <- rep(1.0, length(bnames2))
  init_bp2[which(grepl("lambda", bnames2))] <- rate2
  
  init_fp <- if (frailty %in% c("gamma", "invgaussian")) {
    c(0.5, 0.4)
  } else if (frailty == "posstable") {
    c(0.6, 0.4)
  } else if (frailty == "comppoisson") {
    c(0.5, 1.0, 1.0)
  } else {
    numeric(0)
  }
  
  init_beta <- if (p_cov > 0) rep(0.0, p_cov) else numeric(0)
  init_theta <- c(init_bp1, init_bp2, init_fp, init_beta)
  init_phi <- to_unconstrained(init_theta)
  
  opt_method <- if (!is.null(control$optim.method)) control$optim.method else "BFGS"
  opt_res <- stats::optim(init_phi, neg_loglik_obj, method = opt_method,
                          hessian = TRUE, control = control$optim.control)
  
  mle_phi <- opt_res$par
  mle_theta <- to_constrained(mle_phi)
  names(mle_theta) <- all_param_names
  
  loglik_val <- -opt_res$value
  dev_val <- 2 * opt_res$value
  aic_val <- dev_val + 2 * n_params
  bic_val <- dev_val + n_params * log(n_obs)
  
  # Hessian / Covariance on phi transformed space then Delta method for theta
  vcov_theta <- matrix(NA_real_, nrow = n_params, ncol = n_params)
  colnames(vcov_theta) <- rownames(vcov_theta) <- all_param_names
  std_errs <- rep(NA_real_, n_params)
  names(std_errs) <- all_param_names
  
  tryCatch({
    if (!is.null(opt_res$hessian)) {
      inv_hess <- solve(opt_res$hessian)
      # Numerical gradient d(theta)/d(phi)
      eps <- 1e-5
      grad_mat <- matrix(0, nrow = n_params, ncol = n_params)
      for (k in seq_len(n_params)) {
        phi_plus <- mle_phi; phi_plus[k] <- phi_plus[k] + eps
        phi_minus <- mle_phi; phi_minus[k] <- phi_minus[k] - eps
        grad_mat[, k] <- (to_constrained(phi_plus) - to_constrained(phi_minus)) / (2 * eps)
      }
      vcov_theta <- grad_mat %*% inv_hess %*% t(grad_mat)
      std_errs <- sqrt(pmax(1e-12, diag(vcov_theta)))
    }
  }, error = function(e) NULL)
  
  summary_tab <- data.frame(
    Estimate = mle_theta,
    Std.Error = std_errs,
    z.value = mle_theta / std_errs,
    p.value = 2 * (1 - stats::pnorm(abs(mle_theta / std_errs))),
    Lower.95 = mle_theta - 1.96 * std_errs,
    Upper.95 = mle_theta + 1.96 * std_errs,
    row.names = all_param_names
  )
  
  list(
    method = "mle",
    coefficients = mle_theta,
    vcov = vcov_theta,
    summary.table = summary_tab,
    deviance = dev_val,
    logLik = loglik_val,
    AIC = aic_val,
    BIC = bic_val,
    DIC = NA_real_,
    LPML = NA_real_,
    pD = NA_real_,
    convergence = opt_res$convergence,
    n.params = n_params,
    param.names = all_param_names,
    unpack = unpack_params,
    calc_loglik_indiv = calc_loglik_indiv
  )
}
