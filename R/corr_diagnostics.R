#' MCMC Convergence Diagnostics and Model Selection Criteria
#'
#' @description
#' Functions to calculate Geweke convergence diagnostics, Gelman-Rubin potential
#' scale reduction factors, sample autocorrelations, and model comparison criteria
#' (AIC, BIC, DIC, LPML).
#'
#' @name diagnostics
#' @keywords internal
NULL

corr_calc_geweke <- function(chain, frac1 = 0.1, frac2 = 0.5) {
  chain <- as.numeric(chain)
  n <- length(chain)
  n1 <- floor(n * frac1)
  n2 <- floor(n * frac2)
  
  if (n1 < 2 || n2 < 2) {
    return(list(z = 0.0, p_value = 1.0))
  }
  
  x1 <- chain[1:n1]
  x2 <- chain[(n - n2 + 1):n]
  
  m1 <- mean(x1)
  m2 <- mean(x2)
  
  # Spectral density at zero estimate using AR(1) or sample variance
  v1 <- var(x1) / n1
  v2 <- var(x2) / n2
  
  denom <- sqrt(v1 + v2)
  if (denom <= 1e-12 || is.na(denom)) {
    z <- 0.0
    p_val <- 1.0
  } else {
    z <- (m1 - m2) / denom
    p_val <- 2 * (1 - pnorm(abs(z)))
  }
  
  list(z = z, p_value = p_val)
}

corr_calc_gelman_rubin <- function(chains_list) {
  # chains_list is a list of numeric vectors (one vector per chain for a single parameter)
  m <- length(chains_list)
  if (m < 2) {
    return(1.0)
  }
  
  n <- min(vapply(chains_list, length, integer(1)))
  if (n < 4) {
    return(1.0)
  }
  
  mat <- do.call(cbind, lapply(chains_list, function(ch) ch[1:n]))
  
  chain_means <- colMeans(mat)
  grand_mean <- mean(chain_means)
  
  B <- (n / (m - 1)) * sum((chain_means - grand_mean)^2)
  
  chain_vars <- apply(mat, 2, var)
  W <- mean(chain_vars)
  
  if (W <= 1e-12 || is.na(W)) {
    return(1.0)
  }
  
  var_hat <- ((n - 1) / n) * W + (1 / n) * B
  r_hat <- sqrt(max(1.0, var_hat / W))
  
  if (is.na(r_hat) || is.infinite(r_hat)) 1.0 else r_hat
}

calc_autocorr <- function(chain, lags = c(1, 5, 10, 20, 50, 100)) {
  chain <- as.numeric(chain)
  n <- length(chain)
  chain_centered <- chain - mean(chain)
  c0 <- sum(chain_centered^2)
  
  if (c0 <= 1e-12) {
    res <- rep(0, length(lags))
    names(res) <- paste0("lag_", lags)
    return(res)
  }
  
  res <- numeric(length(lags))
  for (i in seq_along(lags)) {
    k <- lags[i]
    if (k < n) {
      ck <- sum(chain_centered[1:(n - k)] * chain_centered[(k + 1):n])
      res[i] <- ck / c0
    } else {
      res[i] <- 0
    }
  }
  names(res) <- paste0("lag_", lags)
  res
}

corr_calc_model_criteria <- function(loglik_mat, deviance_at_mean, n_obs, p_count) {
  # loglik_mat: matrix of dimension (S iterations x n individuals)
  # deviance_at_mean: D(theta_bar) = -2 * logL(theta_bar)
  
  S <- nrow(loglik_mat)
  n <- ncol(loglik_mat)
  
  # Total log-likelihood for each MCMC draw
  tot_loglik <- rowSums(loglik_mat)
  deviances <- -2 * tot_loglik
  
  mean_deviance <- mean(deviances)
  pD <- max(0, mean_deviance - deviance_at_mean)
  dic <- deviance_at_mean + 2 * pD
  aic <- deviance_at_mean + 2 * p_count
  bic <- deviance_at_mean + p_count * log(n_obs)
  
  # LPML via log-sum-exp harmonic mean
  # CPO_j = S / sum(1 / L_j) -> log(CPO_j) = -log(mean(exp(-loglik_j)))
  lpml_vals <- numeric(n)
  for (j in seq_len(n)) {
    ll_j <- loglik_mat[, j]
    max_ll <- max(ll_j)
    # log(mean(exp(-(ll_j - max_ll)))) - max_ll
    shifted_exp <- exp(-(ll_j - max_ll))
    mean_shifted <- mean(shifted_exp)
    if (mean_shifted > 0 && is.finite(mean_shifted)) {
      lpml_vals[j] <- max_ll - log(mean_shifted)
    } else {
      lpml_vals[j] <- max_ll
    }
  }
  lpml <- sum(lpml_vals)
  
  list(
    deviance_at_mean = deviance_at_mean,
    mean_deviance = mean_deviance,
    pD = pD,
    DIC = dic,
    AIC = aic,
    BIC = bic,
    LPML = lpml
  )
}
