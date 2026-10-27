#' Bayesian Estimation of Shared, Additive, and Multiplicative Frailty Models
#'
#' Fits parametric shared multiplicative, shared additive, univariate multiplicative,
#' univariate additive, and merged frailty survival models using adaptive Markov Chain
#' Monte Carlo (MCMC) sampling.
#'
#' @param formula A formula expression, e.g. \code{Surv(time1, time2, status1, status2) ~ x1 + x2}
#'   for bivariate data, or \code{Surv(time, status) ~ x1 + x2} for univariate data.
#' @param data A data.frame containing all survival variables and covariates.
#' @param baseline Character specifying baseline distribution: \code{"weibull"} (default),
#'   \code{"exponential"}, \code{"genexponential"}, \code{"exppower"}, \code{"lfr"},
#'   \code{"genloglogistic"}, \code{"genweibull"}, \code{"genrayleigh"},
#'   \code{"weightedexp"}, or \code{"extendedweibull"}.
#' @param frailty Character specifying frailty distribution: \code{"gamma"} (default),
#'   \code{"invgauss"}, \code{"posstable"}, \code{"compoundpoisson"},
#'   \code{"compoundnb"}, \code{"genlindley"}, \code{"weightedlindley"},
#'   or \code{"mixture_invgauss"}.
#' @param type Model formulation: \code{"multiplicative_bivariate"} (default),
#'   \code{"additive_bivariate"}, \code{"multiplicative_univariate"},
#'   \code{"additive_univariate"}, or \code{"merged_univariate"}.
#' @param n_iter Number of MCMC iterations saved post burn-in (default 3000).
#' @param burn_in Number of burn-in iterations (default 1000).
#' @param thin Thinning interval (default 1).
#' @param priors List of prior distribution specifications (optional).
#' @param init Vector of initial parameter values (optional).
#' @param adapt_freq Frequency of proposal covariance matrix adaptation (default 50).
#' @param verbose Logical; if \code{TRUE}, displays progress.
#' @param ... Additional arguments.
#'
#' @return An S3 object of class \code{"frailty_bayes"} containing:
#'   \item{chain}{Matrix of posterior draws for all parameters.}
#'   \item{summary}{Data frame of posterior summaries, credible intervals, ESS, MCSE, and Geweke diagnostics.}
#'   \item{diagnostics}{List containing DIC, EAIC, EBIC, CPO, and LPML model comparison criteria.}
#'   \item{log_lik_matrix}{Matrix of individual log-likelihoods across iterations.}
#'   \item{accept_rate}{MCMC acceptance rate.}
#'   \item{data}{Processed data list.}
#'   \item{call}{Matched call.}
#'   \item{baseline}{Baseline distribution name.}
#'   \item{frailty}{Frailty distribution name.}
#'   \item{type}{Model type.}
#'
#' @references
#' Hanagal, D. D. and Dabade, A. D. (2013). Compound negative binomial shared frailty models
#' for bivariate survival data. \emph{Statistics and Probability Letters}, 83(11): 2507-2515.
#' \doi{10.1016/j.spl.2013.07.015}
#'
#' Tyagi, S., Pandey, A., Agiwal, V. and Chesneau, C. (2021). Weighted Lindley multiplicative
#' regression frailty models under random censored data. \emph{Computational and Applied Mathematics},
#' 40: 265. \doi{10.1007/s40314-021-01666-5}
#'
#' Kiprotich, G., Tyagi, S. and Ramos, P. L. (2025). Bayesian survival modeling with mixtures of
#' inverse Gaussian frailties. \emph{Journal of Applied Statistics}. \doi{10.1080/02664763.2025.2560510}
#'
#' @export
#' @examples
#' \donttest{
#' data(kidney_frailty)
#' fit <- frailty_bayes(
#'   Surv(time1, time2, status1, status2) ~ age + sex,
#'   data = kidney_frailty,
#'   baseline = "weibull",
#'   frailty = "gamma",
#'   type = "multiplicative_bivariate",
#'   n_iter = 500,
#'   burn_in = 100
#' )
#' print(fit)
#' summary(fit)
#' }
frailty_bayes <- function(formula, data, baseline = c("weibull", "exponential",
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
                          n_iter = 3000, burn_in = 1000, thin = 1,
                          priors = NULL, init = NULL, adapt_freq = 50,
                          verbose = FALSE, ...) {
  cl <- match.call()
  baseline <- match.arg(baseline)
  frailty <- match.arg(frailty)
  type <- match.arg(type)
  
  # Parse formula and data
  parsed <- parse_frailty_data(formula, data, type)
  n_base <- get_n_base_params(baseline)
  n_frail <- get_n_frail_params(frailty)
  n_cov <- parsed$p_cov
  n_total <- n_base + n_cov + n_frail
  
  param_names <- get_param_names(baseline, frailty, parsed$covar_names)
  
  # Default initial values if not provided
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
  
  # Transformation functions (unconstrained <-> constrained)
  trans_info <- get_param_transforms(baseline, frailty, n_cov)
  to_uncon <- trans_info$to_uncon
  to_con <- trans_info$to_con
  log_jacobian <- trans_info$log_jacobian
  
  init_uncon <- to_uncon(init)
  
  # Define log-posterior function
  log_post_fn <- function(uncon_p) {
    con_p <- to_con(uncon_p)
    if (any(!is.finite(con_p))) return(-1e10)
    
    # Log-likelihood
    ll <- loglik_frailty(con_p, parsed, baseline = baseline,
                         frailty = frailty, type = type, return_indiv = FALSE)
    if (!is.finite(ll) || ll < -1e9) return(-1e10)
    
    # Log-prior
    lp <- eval_log_prior(con_p, n_base, n_cov, n_frail, baseline, frailty, priors)
    
    # Jacobian adjustment
    jac <- log_jacobian(con_p)
    
    total_lp <- ll + lp + jac
    if (!is.finite(total_lp)) return(-1e10)
    total_lp
  }
  
  # Run Adaptive MCMC
  mcmc_res <- adaptive_metropolis(log_post_fn, init_uncon, n_iter = n_iter,
                                  burn_in = burn_in, thin = thin,
                                  adapt_freq = adapt_freq, verbose = verbose)
  
  # Transform chain to constrained space
  chain_con <- t(apply(mcmc_res$chain_uncon, 1, to_con))
  colnames(chain_con) <- param_names
  
  # Compute diagnostics and individual likelihoods
  diag_summary <- compute_mcmc_diagnostics(chain_con, param_names)
  
  # Individual log-likelihood matrix for DIC, CPO, LPML
  log_lik_mat <- matrix(0, nrow = n_iter, ncol = parsed$n)
  for (s in 1:n_iter) {
    log_lik_mat[s, ] <- loglik_frailty(chain_con[s, ], parsed,
                                       baseline = baseline, frailty = frailty,
                                       type = type, return_indiv = TRUE)
  }
  
  # Model comparison criteria
  mean_params <- colMeans(chain_con)
  D_hat <- -2 * loglik_frailty(mean_params, parsed, baseline = baseline,
                               frailty = frailty, type = type, return_indiv = FALSE)
  D_bar <- -2 * mean(rowSums(log_lik_mat))
  p_D <- D_bar - D_hat
  dic_val <- D_bar + p_D
  
  p_count <- n_total
  eaic_val <- D_bar + 2 * p_count
  ebic_val <- D_bar + p_count * log(parsed$n)
  
  # CPO and LPML calculation (Harmonic mean with log-sum-exp)
  # CPO_i = 1 / ( (1/S) * sum_{s=1}^S exp(-log_lik_mat[s, i]) )
  cpo_vec <- numeric(parsed$n)
  for (i in 1:parsed$n) {
    neg_ll <- -log_lik_mat[, i]
    max_neg <- max(neg_ll)
    log_harm_mean <- max_neg + log(mean(exp(neg_ll - max_neg)))
    cpo_vec[i] <- exp(-log_harm_mean)
  }
  lpml_val <- sum(log(pmax(cpo_vec, 1e-300)))
  
  fit_crit <- list(
    DIC = dic_val,
    D_bar = D_bar,
    D_hat = D_hat,
    p_D = p_D,
    EAIC = eaic_val,
    EBIC = ebic_val,
    LPML = lpml_val,
    CPO = cpo_vec,
    LogLik = -0.5 * D_hat,
    n_params = p_count,
    n_obs = parsed$n
  )
  
  out <- list(
    chain = chain_con,
    summary = diag_summary,
    diagnostics = fit_crit,
    log_lik_matrix = log_lik_mat,
    accept_rate = mcmc_res$accept_rate,
    data = parsed,
    call = cl,
    baseline = baseline,
    frailty = frailty,
    type = type,
    param_names = param_names,
    p_cov = n_cov,
    p_base = n_base,
    p_frail = n_frail
  )
  class(out) <- "frailty_bayes"
  out
}
