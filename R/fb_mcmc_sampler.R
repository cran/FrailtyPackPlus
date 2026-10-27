#' Adaptive Metropolis-Hastings MCMC Sampler for Frailty Models
#'
#' Implements an adaptive random-walk Metropolis-Hastings MCMC algorithm
#' targeting an optimal acceptance rate with parameter transformations.
#'
#' @param log_post_fn Function evaluating the log-posterior given unconstrained parameters.
#' @param init_uncon Initial vector of unconstrained parameters.
#' @param n_iter Number of MCMC iterations.
#' @param burn_in Number of burn-in iterations.
#' @param thin Thinning interval.
#' @param adapt_freq Frequency of proposal covariance updates.
#' @param target_accept Target acceptance rate (default 0.234 for d >= 2, 0.44 for d = 1).
#' @param verbose Logical; if \code{TRUE}, prints progress.
#'
#' @return A list containing \code{chain} (unconstrained samples), \code{accept_rate},
#'   \code{scale_factor}, and \code{proposal_cov}.
#' @export
#' @keywords internal
adaptive_metropolis <- function(log_post_fn, init_uncon, n_iter = 5000, burn_in = 1000,
                               thin = 1, adapt_freq = 50, target_accept = NULL,
                               verbose = FALSE) {
  d <- length(init_uncon)
  if (is.null(target_accept)) {
    target_accept <- if (d == 1) 0.44 else 0.234
  }
  
  tot_samples <- burn_in + n_iter * thin
  chain_uncon <- matrix(0, nrow = n_iter, ncol = d)
  log_post_vals <- numeric(n_iter)
  
  curr_uncon <- init_uncon
  curr_log_post <- log_post_fn(curr_uncon)
  
  # Ensure valid initial point
  if (!is.finite(curr_log_post) || curr_log_post < -1e9) {
    # Try slight jitter to find valid point
    for (attempt in 1:100) {
      test_uncon <- init_uncon + stats::rnorm(d, 0, 0.1)
      test_lp <- log_post_fn(test_uncon)
      if (is.finite(test_lp) && test_lp > -1e9) {
        curr_uncon <- test_uncon
        curr_log_post <- test_lp
        break
      }
    }
  }
  
  # Proposal covariance and scaling factor
  prop_cov <- diag(0.01, d)
  log_scale <- log(2.38 / sqrt(d))
  
  # Running mean and covariance for adaptation (Haario et al., 2001)
  running_mean <- curr_uncon
  running_cov <- diag(0.01, d)
  
  n_accepted <- 0
  sample_idx <- 0
  eps_cov <- 1e-6
  
  for (t in 1:tot_samples) {
    scale <- exp(log_scale)
    # Generate proposal
    prop_uncon <- as.vector(curr_uncon + scale * MASS::mvrnorm(1, mu = rep(0, d), Sigma = prop_cov))
    prop_log_post <- log_post_fn(prop_uncon)
    
    # Acceptance probability
    log_alpha <- prop_log_post - curr_log_post
    if (is.na(log_alpha)) log_alpha <- -Inf
    accept_prob <- min(1, exp(pmin(log_alpha, 0)))
    
    if (stats::runif(1) < accept_prob) {
      curr_uncon <- prop_uncon
      curr_log_post <- prop_log_post
      n_accepted <- n_accepted + 1
    }
    
    # Adaptive scale tuning (Robbins-Monro)
    gamma_step <- min(0.1, 1 / (t^0.6))
    log_scale <- log_scale + gamma_step * (accept_prob - target_accept)
    log_scale <- pmin(pmax(log_scale, -10), 10)
    
    # Adaptive covariance update (after initial 100 iterations)
    if (t > 100) {
      diff_m <- curr_uncon - running_mean
      running_mean <- running_mean + (1 / t) * diff_m
      running_cov <- ((t - 1) / t) * running_cov + (1 / t) * (diff_m %*% t(diff_m))
      
      if (t %% adapt_freq == 0) {
        prop_cov <- running_cov + diag(eps_cov, d)
        # Symmetrize and check positive definiteness
        prop_cov <- (prop_cov + t(prop_cov)) / 2
        eig <- eigen(prop_cov, symmetric = TRUE)
        if (any(eig$values <= 0)) {
          prop_cov <- prop_cov + diag(abs(min(eig$values)) + 1e-4, d)
        }
      }
    }
    
    # Save post burn-in and thinned samples
    if (t > burn_in && ((t - burn_in) %% thin == 0)) {
      sample_idx <- sample_idx + 1
      if (sample_idx <= n_iter) {
        chain_uncon[sample_idx, ] <- curr_uncon
        log_post_vals[sample_idx] <- curr_log_post
      }
    }
    
    if (verbose && (t %% 1000 == 0 || t == tot_samples)) {
      utils::setTxtProgressBar(utils::txtProgressBar(min = 0, max = tot_samples, style = 3), t)
    }
  }
  
  accept_rate <- n_accepted / tot_samples
  
  list(chain_uncon = chain_uncon, log_post = log_post_vals,
       accept_rate = accept_rate, prop_cov = prop_cov,
       scale = exp(log_scale))
}
