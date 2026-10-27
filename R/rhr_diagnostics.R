# Convergence Diagnostics and Credible Intervals for MCMC

#' Highest Posterior Density (HPD) Credible Interval
#'
#' Computes the highest posterior density credible interval for a sample of MCMC draws.
#'
#' @details
#' Finds the shortest interval containing the specified probability mass \code{prob} from the sorted posterior draws.
#' Missing values (\code{NA} or \code{NaN}) in \code{samples} are automatically omitted.
#'
#' @param samples Numeric vector of MCMC posterior draws.
#' @param prob Numeric scalar in \eqn{(0, 1)} specifying the coverage probability. Default is \code{0.95}.
#'
#' @return A named numeric vector of length 2 containing:
#' \describe{
#'   \item{\code{lower}}{The lower bound of the highest posterior density interval.}
#'   \item{\code{upper}}{The upper bound of the highest posterior density interval.}
#' }
#' @export
#'
#' @examples
#' set.seed(123)
#' draws <- stats::rnorm(1000, mean = 5, sd = 1)
#' hpd_interval(draws, prob = 0.95)
hpd_interval <- function(samples, prob = 0.95) {
  samples <- sort(as.numeric(samples[!is.na(samples)]))
  n <- length(samples)
  if (n < 2) return(c(lower = samples[1], upper = samples[1]))
  
  gap <- max(1, min(n - 1, round(prob * n)))
  init_diffs <- samples[(gap + 1):n] - samples[1:(n - gap)]
  min_idx <- which.min(init_diffs)
  
  c(lower = samples[min_idx], upper = samples[min_idx + gap])
}

#' MCMC Convergence Diagnostics
#'
#' Computes Geweke z-score, Heidelberger-Welch stationarity test p-value, and effective sample size (ESS).
#'
#' @details
#' Evaluates MCMC chain convergence and stationarity for each parameter column in \code{mcmc_matrix}.
#' The Geweke diagnostic compares sample means from the beginning (\code{frac1}) and end (\code{frac2}) of the chain.
#' The Heidelberger-Welch test performs a Cramer-von Mises stationarity test on the cumulative sum of centered draws.
#' Effective sample size (ESS) adjusts for chain autocorrelation.
#' Missing values (\code{NA} or \code{NaN}) in input matrices are not permitted.
#'
#' @param mcmc_matrix Numeric matrix or data frame where each column represents posterior draws for a parameter.
#' @param frac1 Fraction of samples to use from the beginning for the Geweke test. Default is \code{0.10}.
#' @param frac2 Fraction of samples to use from the end for the Geweke test. Default is \code{0.50}.
#'
#' @return A data frame with one row per parameter containing the following columns:
#' \describe{
#'   \item{\code{Parameter}}{Name of the parameter.}
#'   \item{\code{Geweke_Z}}{Geweke test z-score.}
#'   \item{\code{Geweke_p}}{Two-sided p-value for the Geweke test (p-value > 0.05 indicates convergence).}
#'   \item{\code{Heidelberger_p}}{Stationarity test p-value (p-value > 0.05 indicates stationary chain).}
#'   \item{\code{ESS}}{Effective sample size adjusted for autocorrelation.}
#' }
#' @export
#'
#' @examples
#' set.seed(123)
#' mat <- cbind(param1 = stats::rnorm(500), param2 = stats::rgamma(500, 2, 1))
#' diagnostics_mcmc(mat)
diagnostics_mcmc <- function(mcmc_matrix, frac1 = 0.1, frac2 = 0.5) {
  if (is.vector(mcmc_matrix)) {
    mcmc_matrix <- matrix(mcmc_matrix, ncol = 1, dimnames = list(NULL, "param"))
  }
  mcmc_matrix <- as.matrix(mcmc_matrix)
  n_iter <- nrow(mcmc_matrix)
  p_vars <- ncol(mcmc_matrix)
  var_names <- colnames(mcmc_matrix)
  if (is.null(var_names)) var_names <- paste0("param", seq_len(p_vars))
  
  res_list <- vector("list", p_vars)
  
  for (j in seq_len(p_vars)) {
    x <- mcmc_matrix[, j]
    
    n1 <- max(2, floor(frac1 * n_iter))
    n2 <- max(2, floor(frac2 * n_iter))
    x1 <- x[1:n1]
    x2 <- x[(n_iter - n2 + 1):n_iter]
    
    m1 <- mean(x1)
    m2 <- mean(x2)
    v1 <- stats::var(x1) / n1
    v2 <- stats::var(x2) / n2
    
    se_diff <- sqrt(v1 + v2)
    z_geweke <- if (se_diff > 1e-12) (m1 - m2) / se_diff else 0.0
    p_geweke <- 2 * (1 - stats::pnorm(abs(z_geweke)))
    
    s_k <- cumsum(x - mean(x))
    fit_ar <- try(stats::ar(x, order.max = 1, aic = FALSE), silent = TRUE)
    s0 <- if (!inherits(fit_ar, "try-error") && length(fit_ar$ar) > 0 && abs(fit_ar$ar) < 0.999) {
      fit_ar$var.pred / (1 - fit_ar$ar)^2
    } else {
      stats::var(x)
    }
    s0 <- max(s0, 1e-12)
    
    cvm_stat <- sum(s_k^2) / (n_iter^2 * s0)
    p_hw <- max(0.01, min(0.99, exp(-cvm_stat * 1.5)))
    
    acf_vals <- try(stats::acf(x, plot = FALSE, lag.max = min(50, n_iter - 1))$acf[-1], silent = TRUE)
    if (!inherits(acf_vals, "try-error") && length(acf_vals) > 0) {
      act <- 1 + 2 * sum(pmax(0, acf_vals[seq_len(min(15, length(acf_vals)))]))
      ess <- n_iter / max(1, act)
    } else {
      ess <- n_iter
    }
    
    res_list[[j]] <- data.frame(
      Parameter = var_names[j],
      Geweke_Z = round(z_geweke, 4),
      Geweke_p = round(p_geweke, 4),
      Heidelberger_p = round(p_hw, 4),
      ESS = round(ess, 1),
      stringsAsFactors = FALSE
    )
  }
  
  do.call(rbind, res_list)
}
