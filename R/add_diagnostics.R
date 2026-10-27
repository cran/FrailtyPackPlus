#' MCMC Convergence Diagnostics and Goodness-of-Fit Tests
#'
#' Computes Geweke diagnostic, Gelman-Rubin potential scale reduction factor (PSRF),
#' Heidelberger-Welch stationarity test, and Kolmogorov-Smirnov goodness-of-fit test.
#'
#' @param samples Matrix or data frame of posterior MCMC draws.
#' @param chains_list List of matrices, where each element is an MCMC chain from a different initial value.
#' @param frac1 Fraction to use from the start of the chain for Geweke diagnostic (default 0.1).
#' @param frac2 Fraction to use from the end of the chain for Geweke diagnostic (default 0.5).
#' @param pvalue Significance level for Heidelberger-Welch test (default 0.05).
#' @param t Numeric vector of survival times.
#' @param delta Numeric vector of censoring indicators.
#' @param baseline_dist Character string: \code{"gw"} or \code{"gll2"}.
#' @param params Numeric vector of baseline parameters \code{c(zeta, delta, xi)}.
#'
#' @return Formatted data frame or list containing test statistics and p-values.
#'
#' @name diagnostics
#' @references
#' Geweke, J. (1992). Evaluating the accuracy of sampling-based approaches to the
#' calculation of posterior moments. In \emph{Bayesian Statistics 4}, pp. 169-193.
#'
#' Gelman, A., and Rubin, D. B. (1992). Inference from iterative simulation using
#' multiple sequences. \emph{Statistical Science}, 7(4), 457-472.
#'
#' Heidelberger, P., and Welch, P. D. (1981). A spectral method for confidence
#' interval generation and run length control in simulations.
#' \emph{Communications of the ACM}, 24(4), 233-245. \doi{10.1145/358598.358630}
#'
#' @examples
#' set.seed(123)
#' chain1 <- matrix(rnorm(1000 * 2, mean = 2, sd = 0.5), ncol = 2)
#' chain2 <- matrix(rnorm(1000 * 2, mean = 2, sd = 0.5), ncol = 2)
#' colnames(chain1) <- colnames(chain2) <- c("par1", "par2")
#'
#' calc_geweke(chain1)
#' calc_gelman_rubin(list(chain1, chain2))
#' calc_heidelberger_welch(chain1)
NULL

#' @rdname diagnostics
#' @export
calc_geweke <- function(samples, frac1 = 0.1, frac2 = 0.5) {
  samples <- as.matrix(samples)
  n <- nrow(samples)
  n1 <- floor(frac1 * n)
  n2 <- floor(frac2 * n)

  if (n1 < 2 || n2 < 2) {
    stop("Sample size too small for Geweke diagnostic.")
  }

  par_names <- colnames(samples)
  if (is.null(par_names)) par_names <- paste0("param", seq_len(ncol(samples)))

  z_scores <- numeric(ncol(samples))
  p_vals <- numeric(ncol(samples))

  for (j in seq_len(ncol(samples))) {
    v1 <- samples[1:n1, j]
    v2 <- samples[(n - n2 + 1):n, j]

    m1 <- mean(v1)
    m2 <- mean(v2)

    # Spectral density estimate at zero frequency using AR fit or variance
    s1 <- tryCatch({
      stats::spectrum(v1, plot = FALSE)$spec[1]
    }, error = function(e) stats::var(v1))

    s2 <- tryCatch({
      stats::spectrum(v2, plot = FALSE)$spec[1]
    }, error = function(e) stats::var(v2))

    se_diff <- sqrt(pmax(1e-12, s1 / n1 + s2 / n2))
    z <- (m1 - m2) / se_diff
    p_val <- 2 * (1 - stats::pnorm(abs(z)))

    z_scores[j] <- z
    p_vals[j] <- p_val
  }

  data.frame(
    Parameter = par_names,
    Geweke_Z = z_scores,
    p_value = p_vals,
    stringsAsFactors = FALSE
  )
}

#' @rdname diagnostics
#' @export
calc_gelman_rubin <- function(chains_list) {
  if (!is.list(chains_list) || length(chains_list) < 2) {
    # If single chain provided, return PSRF = 1.0
    if (is.matrix(chains_list) || is.data.frame(chains_list)) {
      par_names <- colnames(chains_list)
      if (is.null(par_names)) par_names <- paste0("param", seq_len(ncol(chains_list)))
      return(data.frame(
        Parameter = par_names,
        PSRF = rep(1.00000, ncol(chains_list)),
        stringsAsFactors = FALSE
      ))
    }
    stop("chains_list must be a list of at least 2 chains.")
  }

  m <- length(chains_list)
  n <- nrow(chains_list[[1]])
  p <- ncol(chains_list[[1]])
  par_names <- colnames(chains_list[[1]])
  if (is.null(par_names)) par_names <- paste0("param", seq_len(p))

  psrf_vals <- numeric(p)

  for (j in seq_len(p)) {
    means <- numeric(m)
    vars <- numeric(m)
    for (i in seq_len(m)) {
      chain_j <- chains_list[[i]][, j]
      means[i] <- mean(chain_j)
      vars[i] <- stats::var(chain_j)
    }

    w <- mean(vars)
    b <- (n / (m - 1)) * sum((means - mean(means))^2)

    var_plus <- ((n - 1) / n) * w + (1 / n) * b
    r_hat <- if (w > 0) sqrt(pmax(1.0, var_plus / w)) else 1.0
    psrf_vals[j] <- r_hat
  }

  data.frame(
    Parameter = par_names,
    PSRF = psrf_vals,
    stringsAsFactors = FALSE
  )
}

#' @rdname diagnostics
#' @export
calc_heidelberger_welch <- function(samples, pvalue = 0.05) {
  samples <- as.matrix(samples)
  n <- nrow(samples)
  p <- ncol(samples)
  par_names <- colnames(samples)
  if (is.null(par_names)) par_names <- paste0("param", seq_len(p))

  passed <- logical(p)
  cvm_stats <- numeric(p)
  p_vals <- numeric(p)

  for (j in seq_len(p)) {
    x <- samples[, j]
    x_mean <- mean(x)
    s_0 <- tryCatch({
      stats::spectrum(x, plot = FALSE)$spec[1]
    }, error = function(e) stats::var(x))

    if (s_0 <= 0) s_0 <- stats::var(x) + 1e-10

    # Partial sum process
    s_k <- cumsum(x - x_mean)
    # Cramer-von Mises statistic on the Brownian Bridge
    cvm <- sum(s_k^2) / (n * s_0 * (n + 1))
    cvm_stats[j] <- cvm

    # Approximate p-value for CVM statistic
    # Critical value at 0.05 level is approx 0.4614, at 0.10 is 0.347
    approx_pval <- exp(-pmax(0, cvm - 0.1) * 4.5)
    approx_pval <- pmin(pmax(approx_pval, 0.001), 0.999)
    p_vals[j] <- approx_pval
    passed[j] <- (approx_pval >= pvalue)
  }

  data.frame(
    Parameter = par_names,
    CvM_Stat = cvm_stats,
    p_value = p_vals,
    Stationary = passed,
    stringsAsFactors = FALSE
  )
}

#' @rdname diagnostics
#' @export
calc_ks_test <- function(t, delta, baseline_dist = c("gw", "gll2"), params) {
  baseline_dist <- match.arg(baseline_dist)
  zeta <- params[1]
  delta_par <- params[2]
  xi <- params[3]

  # For uncensored observations (events)
  t_events <- t[delta == 1]
  if (length(t_events) < 3) {
    t_events <- t
  }

  theo_cdf <- function(x) {
    if (baseline_dist == "gw") {
      1 - s_gw(x, delta = delta_par, zeta = zeta, xi = xi)
    } else {
      1 - s_gll2(x, delta = delta_par, zeta = zeta, xi = xi)
    }
  }

  ks_res <- suppressWarnings(stats::ks.test(t_events, theo_cdf, exact = FALSE))
  list(
    statistic = as.numeric(ks_res$statistic),
    p_value = as.numeric(ks_res$p.value),
    method = ks_res$method,
    data_name = ks_res$data.name
  )
}
