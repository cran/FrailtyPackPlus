#' MCMC Estimation Engine for Correlated Frailty Models
#'
#' @description
#' Adaptive Metropolis-Hastings within Gibbs sampler for Bayesian estimation
#' of correlated frailty models for bivariate survival data.
#'
#' @name mcmc_engine
#' @keywords internal
NULL

fit_corrFrailty_mcmc <- function(t1, d1, t2, d2, X,
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
  
  # Helper to unpack vector into components
  unpack_params <- function(theta) {
    nb1 <- length(bnames1)
    nb2 <- length(bnames2)
    nf <- length(fnames)
    
    idx <- 1
    bp1 <- theta[idx:(idx + nb1 - 1)]; idx <- idx + nb1
    bp2 <- theta[idx:(idx + nb2 - 1)]; idx <- idx + nb2
    
    fp <- if (nf > 0) {
      res <- theta[idx:(idx + nf - 1)]
      idx <- idx + nf
      res
    } else {
      numeric(0)
    }
    
    beta <- if (p_cov > 0) theta[idx:(idx + p_cov - 1)] else numeric(0)
    
    list(base1 = bp1, base2 = bp2, frailty = fp, beta = beta)
  }
  
  # Log-prior
  log_prior_fn <- function(unpacked) {
    # Baseline 1 & 2: gamma(1, 0.0001) / exponential diffuse prior
    lp <- 0
    for (val in c(unpacked$base1, unpacked$base2)) {
      if (val <= 0) return(-Inf)
      lp <- lp + stats::dgamma(val, shape = 1, rate = 0.001, log = TRUE)
    }
    
    # Frailty parameters
    if (frailty %in% c("gamma", "invgaussian")) {
      sig <- unpacked$frailty[1]
      rho <- unpacked$frailty[2]
      if (sig <= 0 || rho <= 0 || rho >= 1) return(-Inf)
      lp <- lp + stats::dgamma(sig, shape = 1, rate = 0.001, log = TRUE)
      # Uniform(0, 1) for rho -> log(1) = 0
    } else if (frailty == "posstable") {
      theta <- unpacked$frailty[1]
      rho <- unpacked$frailty[2]
      if (theta <= 0 || theta >= 1 || rho <= 0 || rho >= 1) return(-Inf)
      # Uniform(0, 1) for theta and rho
    } else if (frailty == "comppoisson") {
      eta <- unpacked$frailty[1]
      nu <- unpacked$frailty[2]
      theta <- unpacked$frailty[3]
      if (eta <= 0 || nu <= 0 || theta <= 0) return(-Inf)
      lp <- lp + stats::dgamma(eta, shape = 1, rate = 0.001, log = TRUE) +
                 stats::dgamma(nu, shape = 1, rate = 0.001, log = TRUE) +
                 stats::dgamma(theta, shape = 1, rate = 0.001, log = TRUE)
    }
    
    # Covariates: Normal(0, 1000)
    if (p_cov > 0) {
      for (b in unpacked$beta) {
        lp <- lp + stats::dnorm(b, mean = 0, sd = 31.62, log = TRUE)
      }
    }
    
    lp
  }
  
  # Log-likelihood
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
  
  log_post_fn <- function(theta) {
    unp <- unpack_params(theta)
    lprior <- log_prior_fn(unp)
    if (is.infinite(lprior) && lprior < 0) return(-Inf)
    
    ll_vec <- calc_loglik_indiv(unp)
    if (any(is.nan(ll_vec)) || any(is.na(ll_vec))) return(-Inf)
    
    lprior + sum(ll_vec)
  }
  
  # Parameter transformation for unconstrained random-walk
  to_unconstrained <- function(theta) {
    unp <- unpack_params(theta)
    nb1 <- length(bnames1); nb2 <- length(bnames2); nf <- length(fnames)
    
    phi <- numeric(n_params)
    idx <- 1
    
    # Baseline: log
    phi[idx:(idx + nb1 - 1)] <- log(pmax(1e-6, unp$base1)); idx <- idx + nb1
    phi[idx:(idx + nb2 - 1)] <- log(pmax(1e-6, unp$base2)); idx <- idx + nb2
    
    # Frailty
    if (frailty %in% c("gamma", "invgaussian")) {
      phi[idx] <- log(pmax(1e-6, unp$frailty[1])) # sigma
      phi[idx + 1] <- log(unp$frailty[2] / (1 - unp$frailty[2])) # logit(rho)
      idx <- idx + 2
    } else if (frailty == "posstable") {
      phi[idx] <- log(unp$frailty[1] / (1 - unp$frailty[1])) # logit(theta)
      phi[idx + 1] <- log(unp$frailty[2] / (1 - unp$frailty[2])) # logit(rho)
      idx <- idx + 2
    } else if (frailty == "comppoisson") {
      phi[idx:(idx + 2)] <- log(pmax(1e-6, unp$frailty)) # log(eta, nu, theta)
      idx <- idx + 3
    }
    
    # Beta: identity
    if (p_cov > 0) {
      phi[idx:(idx + p_cov - 1)] <- unp$beta
    }
    
    phi
  }
  
  to_constrained <- function(phi) {
    nb1 <- length(bnames1); nb2 <- length(bnames2); nf <- length(fnames)
    theta <- numeric(n_params)
    idx <- 1
    
    # Baseline: exp
    theta[idx:(idx + nb1 - 1)] <- exp(phi[idx:(idx + nb1 - 1)]); idx <- idx + nb1
    theta[idx:(idx + nb2 - 1)] <- exp(phi[idx:(idx + nb2 - 1)]); idx <- idx + nb2
    
    # Frailty
    if (frailty %in% c("gamma", "invgaussian")) {
      theta[idx] <- exp(phi[idx]) # sigma
      theta[idx + 1] <- 1 / (1 + exp(-phi[idx + 1])) # rho in (0, 1)
      idx <- idx + 2
    } else if (frailty == "posstable") {
      theta[idx] <- 1 / (1 + exp(-phi[idx])) # theta in (0, 1)
      theta[idx + 1] <- 1 / (1 + exp(-phi[idx + 1])) # rho in (0, 1)
      idx <- idx + 2
    } else if (frailty == "comppoisson") {
      theta[idx:(idx + 2)] <- exp(phi[idx:(idx + 2)]) # eta, nu, theta
      idx <- idx + 3
    }
    
    # Beta
    if (p_cov > 0) {
      theta[idx:(idx + p_cov - 1)] <- phi[idx:(idx + p_cov - 1)]
    }
    
    theta
  }
  
  # Log Jacobian for transformed parameter space
  log_jacobian <- function(phi) {
    nb1 <- length(bnames1); nb2 <- length(bnames2); nf <- length(fnames)
    lj <- 0
    idx <- 1
    
    # Baseline
    lj <- lj + sum(phi[idx:(idx + nb1 - 1)]); idx <- idx + nb1
    lj <- lj + sum(phi[idx:(idx + nb2 - 1)]); idx <- idx + nb2
    
    # Frailty
    if (frailty %in% c("gamma", "invgaussian")) {
      lj <- lj + phi[idx] # log(sigma)
      p_rho <- 1 / (1 + exp(-phi[idx + 1]))
      lj <- lj + log(pmax(1e-12, p_rho * (1 - p_rho))) # logit jacobian
      idx <- idx + 2
    } else if (frailty == "posstable") {
      p_th <- 1 / (1 + exp(-phi[idx]))
      p_rho <- 1 / (1 + exp(-phi[idx + 1]))
      lj <- lj + log(pmax(1e-12, p_th * (1 - p_th))) + log(pmax(1e-12, p_rho * (1 - p_rho)))
      idx <- idx + 2
    } else if (frailty == "comppoisson") {
      lj <- lj + sum(phi[idx:(idx + 2)])
      idx <- idx + 3
    }
    
    lj
  }
  
  # Initial values generator
  get_initial_values <- function(chain_num) {
    if (!is.null(control$init) && is.list(control$init) && length(control$init) >= chain_num) {
      init_vals <- control$init[[chain_num]]
      if (length(init_vals) == n_params) return(init_vals)
    }
    
    # Reasonable default heuristics based on empirical hazard
    rate_est1 <- max(0.01, sum(d1) / sum(t1))
    rate_est2 <- max(0.01, sum(d2) / sum(t2))
    
    init_bp1 <- rep(1.0, length(bnames1))
    init_bp1[which(grepl("lambda", bnames1))] <- rate_est1
    
    init_bp2 <- rep(1.0, length(bnames2))
    init_bp2[which(grepl("lambda", bnames2))] <- rate_est2
    
    init_fp <- if (frailty %in% c("gamma", "invgaussian")) {
      c(0.5 + 0.1 * (chain_num - 1), 0.4 + 0.1 * (chain_num - 1))
    } else if (frailty == "posstable") {
      c(0.6 + 0.1 * (chain_num - 1), 0.4 + 0.1 * (chain_num - 1))
    } else if (frailty == "comppoisson") {
      c(0.5 + 0.1 * (chain_num - 1), 1.0, 1.0)
    } else {
      numeric(0)
    }
    
    init_beta <- if (p_cov > 0) rep(0.0, p_cov) else numeric(0)
    
    # Slight perturbation for chain > 1
    if (chain_num > 1) {
      init_bp1 <- init_bp1 * stats::runif(length(init_bp1), 0.8, 1.2)
      init_bp2 <- init_bp2 * stats::runif(length(init_bp2), 0.8, 1.2)
    }
    
    c(init_bp1, init_bp2, init_fp, init_beta)
  }
  
  if (!is.null(control$seed)) {
    set.seed(control$seed)
  }
  
  n_chains <- control$n.chains
  n_iter <- control$n.iter
  n_burnin <- control$n.burnin
  n_thin <- control$n.thin
  adapt_iter <- control$adapt.iter
  
  n_saved <- floor((n_iter - n_burnin) / n_thin)
  if (n_saved < 1) stop("n.iter must be strictly greater than n.burnin")
  
  chains_list <- vector("list", n_chains)
  accept_rates <- matrix(0, nrow = n_chains, ncol = n_params)
  colnames(accept_rates) <- all_param_names
  
  loglik_samples_list <- vector("list", n_chains)
  
  for (chain_i in seq_len(n_chains)) {
    cur_theta <- get_initial_values(chain_i)
    cur_phi <- to_unconstrained(cur_theta)
    
    cur_logpost <- log_post_fn(cur_theta) + log_jacobian(cur_phi)
    
    # Proposal standard deviations for unconstrained space
    step_sd <- if (!is.null(control$tune) && length(control$tune) == n_params) {
      control$tune
    } else {
      rep(0.1, n_params)
    }
    
    chain_samples <- matrix(0, nrow = n_saved, ncol = n_params)
    colnames(chain_samples) <- all_param_names
    
    loglik_samples <- matrix(0, nrow = n_saved, ncol = n_obs)
    
    accept_counts <- integer(n_params)
    total_proposals <- integer(n_params)
    save_idx <- 1
    
    for (iter in seq_len(n_iter)) {
      # Component-wise Metropolis-Hastings (Gibbs-style update)
      for (p_idx in seq_len(n_params)) {
        prop_phi <- cur_phi
        prop_phi[p_idx] <- stats::rnorm(1, mean = cur_phi[p_idx], sd = step_sd[p_idx])
        
        prop_theta <- to_constrained(prop_phi)
        prop_logpost <- log_post_fn(prop_theta) + log_jacobian(prop_phi)
        
        log_alpha <- prop_logpost - cur_logpost
        
        total_proposals[p_idx] <- total_proposals[p_idx] + 1L
        if (is.finite(log_alpha) && (log_alpha >= 0 || stats::runif(1) < exp(log_alpha))) {
          cur_phi <- prop_phi
          cur_theta <- prop_theta
          cur_logpost <- prop_logpost
          accept_counts[p_idx] <- accept_counts[p_idx] + 1L
        }
      }
      
      # Adaptive tuning during burn-in
      if (iter <= adapt_iter && (iter %% 50 == 0)) {
        for (p_idx in seq_len(n_params)) {
          acc_rate_p <- accept_counts[p_idx] / total_proposals[p_idx]
          if (acc_rate_p < 0.2) {
            step_sd[p_idx] <- max(0.01, step_sd[p_idx] * 0.8)
          } else if (acc_rate_p > 0.45) {
            step_sd[p_idx] <- min(2.0, step_sd[p_idx] * 1.2)
          }
        }
      }
      
      # Save sample after burn-in with thinning
      if (iter > n_burnin && ((iter - n_burnin) %% n_thin == 0)) {
        chain_samples[save_idx, ] <- cur_theta
        unp_cur <- unpack_params(cur_theta)
        loglik_samples[save_idx, ] <- calc_loglik_indiv(unp_cur)
        save_idx <- save_idx + 1
      }
    }
    
    chains_list[[chain_i]] <- chain_samples
    accept_rates[chain_i, ] <- accept_counts / total_proposals
    loglik_samples_list[[chain_i]] <- loglik_samples
  }
  
  # Combine chains for posterior inference
  all_samples_mat <- do.call(rbind, chains_list)
  all_loglik_mat <- do.call(rbind, loglik_samples_list)
  
  post_mean <- colMeans(all_samples_mat)
  post_sd <- apply(all_samples_mat, 2, stats::sd)
  post_median <- apply(all_samples_mat, 2, stats::median)
  post_q025 <- apply(all_samples_mat, 2, function(x) stats::quantile(x, 0.025))
  post_q975 <- apply(all_samples_mat, 2, function(x) stats::quantile(x, 0.975))
  
  # Effective sample size and Monte Carlo Standard Error
  mc_error <- post_sd / sqrt(nrow(all_samples_mat))
  
  # Diagnostics
  geweke_res <- lapply(seq_len(n_params), function(j) {
    corr_calc_geweke(all_samples_mat[, j])
  })
  geweke_z <- vapply(geweke_res, function(x) x$z, numeric(1))
  geweke_p <- vapply(geweke_res, function(x) x$p_value, numeric(1))
  
  gelman_rhat <- if (n_chains >= 2) {
    vapply(seq_len(n_params), function(j) {
      ch_j <- lapply(chains_list, function(mat) mat[, j])
      corr_calc_gelman_rubin(ch_j)
    }, numeric(1))
  } else {
    rep(1.0, n_params)
  }
  
  autocorr_mat <- t(vapply(seq_len(n_params), function(j) {
    calc_autocorr(all_samples_mat[, j])
  }, numeric(6)))
  rownames(autocorr_mat) <- all_param_names
  
  # Deviance at posterior mean
  unp_mean <- unpack_params(post_mean)
  loglik_at_mean <- sum(calc_loglik_indiv(unp_mean))
  deviance_at_mean <- -2 * loglik_at_mean
  
  crit <- corr_calc_model_criteria(all_loglik_mat, deviance_at_mean, n_obs, n_params)
  
  # Summary table
  summary_tab <- data.frame(
    Estimate = post_mean,
    Std.Error = post_sd,
    MC.Error = mc_error,
    Median = post_median,
    Lower.95 = post_q025,
    Upper.95 = post_q975,
    Geweke.Z = geweke_z,
    Geweke.p = geweke_p,
    R.hat = gelman_rhat,
    row.names = all_param_names
  )
  
  list(
    method = "mcmc",
    coefficients = post_mean,
    summary.table = summary_tab,
    chains = chains_list,
    all.samples = all_samples_mat,
    acceptance.rates = accept_rates,
    autocorrelation = autocorr_mat,
    deviance = deviance_at_mean,
    logLik = loglik_at_mean,
    pD = crit$pD,
    DIC = crit$DIC,
    AIC = crit$AIC,
    BIC = crit$BIC,
    LPML = crit$LPML,
    n.params = n_params,
    param.names = all_param_names,
    unpack = unpack_params,
    calc_loglik_indiv = calc_loglik_indiv
  )
}
