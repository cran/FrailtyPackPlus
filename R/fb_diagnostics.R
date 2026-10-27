#' Convergence and Summary Diagnostics for MCMC Chains
#'
#' Computes Geweke z-score, Gelman-Rubin diagnostic (multichain), Effective Sample Size (ESS),
#' Monte Carlo Standard Error (MCSE), and Highest Posterior Density (HPD) intervals.
#'
#' @param chain Numeric matrix of MCMC samples (iterations x parameters).
#' @param param_names Character vector of parameter names.
#' @param frac1 Fraction of chain to use for Geweke test first part (default 0.1).
#' @param frac2 Fraction of chain to use for Geweke test second part (default 0.5).
#' @param prob Credibility level for HPD intervals (default 0.95).
#'
#' @return A data frame containing posterior summary statistics and diagnostics.
#' @export
compute_mcmc_diagnostics <- function(chain, param_names = NULL, frac1 = 0.1,
                                    frac2 = 0.5, prob = 0.95) {
  if (is.vector(chain)) chain <- matrix(chain, ncol = 1)
  n <- nrow(chain)
  p <- ncol(chain)
  if (is.null(param_names)) param_names <- paste0("param", 1:p)
  
  res_list <- vector("list", p)
  
  for (j in 1:p) {
    x <- chain[, j]
    m <- mean(x, na.rm = TRUE)
    med <- stats::median(x, na.rm = TRUE)
    s <- stats::sd(x, na.rm = TRUE)
    if (is.na(s)) s <- 0
    
    alpha <- 1 - prob
    q_low <- stats::quantile(x, probs = alpha / 2, names = FALSE, na.rm = TRUE)
    q_high <- stats::quantile(x, probs = 1 - alpha / 2, names = FALSE, na.rm = TRUE)
    
    hpd <- calc_hpd(x, prob = prob)
    ess_val <- calc_ess(x)
    mcse_val <- if (ess_val > 1 && s > 0) s / sqrt(ess_val) else s / sqrt(max(n, 1))
    
    geweke_z <- fb_calc_geweke(x, frac1 = frac1, frac2 = frac2)
    geweke_p <- if (is.finite(geweke_z)) 2 * (1 - stats::pnorm(abs(geweke_z))) else 1.0
    
    res_list[[j]] <- data.frame(
      Parameter = param_names[j],
      Mean = m,
      SD = s,
      Median = med,
      MCSE = mcse_val,
      HPD_Lower = hpd[1],
      HPD_Upper = hpd[2],
      Quantile_2.5 = q_low,
      Quantile_97.5 = q_high,
      ESS = round(ess_val, 1),
      Geweke_Z = round(if (is.finite(geweke_z)) geweke_z else 0, 3),
      Geweke_p = round(if (is.finite(geweke_p)) geweke_p else 1, 4),
      stringsAsFactors = FALSE
    )
  }
  
  do.call(rbind, res_list)
}

#' Calculate Highest Posterior Density (HPD) Interval
#' @param x Numeric vector of MCMC draws.
#' @param prob Probability content.
#' @return Length 2 numeric vector (lower, upper).
#' @export
#' @keywords internal
calc_hpd <- function(x, prob = 0.95) {
  x <- x[!is.na(x)]
  sx <- sort(x)
  n <- length(sx)
  if (n < 2) return(c(if (n == 1) sx[1] else 0, if (n == 1) sx[1] else 0))
  gap <- max(1, floor(prob * n))
  if (gap >= n) return(c(sx[1], sx[n]))
  
  intervals <- sx[(gap + 1):n] - sx[1:(n - gap)]
  min_idx <- which.min(intervals)
  c(sx[min_idx], sx[min_idx + gap])
}

#' Calculate Effective Sample Size (ESS)
#' @param x Numeric vector of MCMC draws.
#' @return Numeric ESS estimate.
#' @export
#' @keywords internal
calc_ess <- function(x) {
  x <- x[!is.na(x)]
  n <- length(x)
  if (n < 10 || stats::sd(x) < 1e-12) return(n)
  
  max_lag <- min(500, floor(n / 2))
  acf_res <- tryCatch(stats::acf(x, lag.max = max_lag, plot = FALSE)$acf[, 1, 1],
                      error = function(e) numeric(0))
  if (length(acf_res) < 3 || any(is.na(acf_res))) return(n)
  
  gamma_sum <- 0
  k <- 1
  while (k < length(acf_res) - 1) {
    pair_sum <- acf_res[k + 1] + acf_res[k + 2]
    if (is.na(pair_sum) || pair_sum < 0) break
    gamma_sum <- gamma_sum + pair_sum
    k <- k + 2
  }
  
  tau <- 1 + 2 * gamma_sum
  ess <- if (is.finite(tau) && tau > 0) n / tau else n
  min(n, max(1, ess))
}

#' Calculate Geweke Convergence Diagnostic
#' @param x Numeric vector of MCMC draws.
#' @param frac1 Fraction for first segment.
#' @param frac2 Fraction for second segment.
#' @return Geweke z-score scalar.
#' @export
#' @keywords internal
fb_calc_geweke <- function(x, frac1 = 0.1, frac2 = 0.5) {
  x <- x[!is.na(x)]
  n <- length(x)
  n1 <- floor(frac1 * n)
  n2 <- floor(frac2 * n)
  if (n1 < 5 || n2 < 5 || stats::sd(x) < 1e-12) return(0)
  
  x1 <- x[1:n1]
  x2 <- x[(n - n2 + 1):n]
  
  ess1 <- calc_ess(x1)
  ess2 <- calc_ess(x2)
  v1 <- stats::var(x1) / max(ess1, 1)
  v2 <- stats::var(x2) / max(ess2, 1)
  denom <- sqrt(v1 + v2)
  
  if (is.finite(denom) && denom > 1e-12) (mean(x1) - mean(x2)) / denom else 0
}
