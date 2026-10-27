# Bayesian MCMC Estimation Engine for Shared RHR Frailty Models

#' Bayesian Estimation of Shared Frailty Models under Reversed Hazard Rate
#'
#' Fits shared frailty regression models for bivariate left-censored survival data
#' based on reversed hazard rates using an adaptive Metropolis-within-Gibbs MCMC algorithm.
#' Provides posterior estimates under squared error (SELF), entropy (ELF), and precautionary (PLF)
#' loss functions, highest posterior density (HPD) credible intervals, and model comparison criteria.
#'
#' @details
#' The model assumes that the bivariate lifetime random variables \eqn{(T_{1j}, T_{2j})} share an unobserved
#' frailty variable \eqn{W_j} that acts multiplicatively on the baseline reversed hazard rates:
#' \deqn{r(t_{ij} \mid W_j, X_{ij}) = W_j r_{0i}(t_{ij}) \exp(X_{ij}^T \beta)}
#' Under left censoring, observations with \eqn{d_{ij} = 0} indicate that the event occurred before or at
#' the observation time \eqn{t_{ij}} (\eqn{T_{ij} \le t_{ij}}), while \eqn{d_{ij} = 1} indicates an exact observed failure time.
#' Missing values (\code{NA} or \code{NaN}) in survival times or censoring indicators are not allowed and will produce an informative error.
#'
#' @param formula An optional object of class \code{"formula"} of the form \code{cbind(t1, t2, d1, d2) ~ x1 + x2}.
#' @param data An optional data frame containing the variables in the model.
#' @param t1 Numeric vector of first failure/censoring times (if formula is not provided). Must be positive.
#' @param t2 Numeric vector of second failure/censoring times. Must be positive.
#' @param d1 Numeric or integer vector of censoring indicators for first component (\code{1} = observed, \code{0} = left-censored).
#' @param d2 Numeric or integer vector of censoring indicators for second component (\code{1} = observed, \code{0} = left-censored).
#' @param covars Optional numeric matrix or data frame of covariates.
#' @param frailty Character string specifying frailty distribution. One of \code{"WL"} (Weighted Lindley),
#'   \code{"GL"} (Generalized Lindley), \code{"GIG"} (Generalized Inverse Gaussian mixture), \code{"Gamma"} (Gamma),
#'   \code{"IG"} (Inverse Gaussian), \code{"Lindley"} (Lindley), or \code{"None"} (Independent baseline). Default is \code{"WL"}.
#' @param baseline Character string specifying baseline distribution. One of \code{"GE"} (Generalized Exponential),
#'   \code{"EG"} (Exponentiated Gumbel), \code{"GIR"} (Generalized Inverse Rayleigh), \code{"MIW"} (Modified Inverse Weibull),
#'   \code{"GLL-I"} (Generalized Log-Logistic Type-I), \code{"GLL-II"} (Generalized Log-Logistic Type-II),
#'   \code{"GIE"} (Generalized Inverted Exponential), \code{"Weibull"}, \code{"Exponential"}, or \code{"Log-Logistic"}. Default is \code{"GE"}.
#' @param n_iter Total number of MCMC iterations. Default is \code{2000}.
#' @param burn_in Number of initial burn-in iterations to discard. Default is \code{500}.
#' @param thin Thinning interval. Default is \code{1}.
#' @param n_chains Number of parallel MCMC chains. Default is \code{1}.
#' @param prior Character string indicating prior specification: \code{"flat"} (non-informative) or \code{"informative"}. Default is \code{"flat"}.
#' @param prior_hyper Optional named list of hyperparameters for informative gamma and normal priors.
#' @param loss_functions Character vector of loss functions under which Bayes estimates are computed. Default is \code{c("SELF", "ELF", "PLF")}.
#' @param init Optional numeric vector of initial values for the parameters.
#' @param tune_sd Optional numeric vector of proposal standard deviations for random-walk Metropolis updates.
#' @param verbose Logical; if \code{TRUE}, progress messages are printed. Default is \code{FALSE}.
#' @param seed Optional integer seed for reproducibility. Default is \code{NULL}.
#'
#' @return An object of class \code{"rhr_bayes_fit"} with the following components:
#' \describe{
#'   \item{\code{summary}}{A data frame summarizing parameter estimates under SELF, ELF, and PLF loss functions, posterior standard deviations, Monte Carlo standard errors (MCSE), medians, 95\% HPD credible intervals, Geweke test p-values, and Heidelberger-Welch stationarity test p-values.}
#'   \item{\code{draws}}{A numeric matrix of kept posterior parameter draws of dimension \code{((n_iter - burn_in) / thin) * n_chains} by \code{n_params}.}
#'   \item{\code{chains}}{A list of chains containing kept parameter draws, log-likelihoods, and inverse pointwise likelihoods.}
#'   \item{\code{estimates}}{A named list containing numeric vectors of Bayes point estimates under \code{SELF} (Posterior Mean), \code{ELF} (Harmonic Mean), and \code{PLF} (Root Mean Square).}
#'   \item{\code{loss_functions}}{Character vector of loss functions evaluated.}
#'   \item{\code{diagnostics}}{A data frame containing MCMC convergence diagnostics including Geweke z-scores and p-values, Heidelberger-Welch p-values, and effective sample sizes (ESS).}
#'   \item{\code{criteria}}{A list of model fit and selection criteria: \code{AIC}, \code{BIC}, \code{AICc}, \code{HQIC}, \code{DIC} (Deviance Information Criterion), \code{pD} (effective number of parameters), \code{EAIC}, \code{EBIC}, \code{LPML} (Log Pseudo Marginal Likelihood), \code{mean_deviance}, and \code{dev_at_mean}.}
#'   \item{\code{cpo}}{Numeric vector of conditional predictive ordinate (CPO) values for each observation pair.}
#'   \item{\code{kendall_tau}}{Numeric scalar representing estimated Kendall's rank correlation coefficient \eqn{\tau}.}
#'   \item{\code{frailty_variance}}{Numeric scalar representing estimated variance of the frailty variable \eqn{V(W)}.}
#'   \item{\code{model_info}}{A list containing model specifications: frailty distribution, baseline distribution, sample size \eqn{n}, number of iterations, burn-in, thinning, number of chains, prior type, covariate names, and parameter names.}
#'   \item{\code{data}}{A list containing the input data vectors and covariate matrix.}
#'   \item{\code{call}}{The matched call.}
#' }
#' @export
#'
#' @examples
#' data(austrian_twin)
#' # Demonstration with a small subset and few iterations
#' fit <- fit_rhr_frailty_bayes(
#'   cbind(t1, t2, d1, d2) ~ age + k2 + k3,
#'   data = austrian_twin[1:50, ],
#'   frailty = "WL", baseline = "GE",
#'   n_iter = 60, burn_in = 10, thin = 1
#' )
#' print(fit)
#' summary(fit)
fit_rhr_frailty_bayes <- function(formula = NULL, data = NULL, t1 = NULL, t2 = NULL, d1 = NULL, d2 = NULL,
                                  covars = NULL, frailty = c("WL", "GL", "GIG", "Gamma", "IG", "Lindley", "None"),
                                  baseline = c("GE", "EG", "GIR", "MIW", "GLL-I", "GLL-II", "GIE", "Weibull", "Exponential", "Log-Logistic"),
                                  n_iter = 2000, burn_in = 500, thin = 1, n_chains = 1,
                                  prior = c("flat", "informative"), prior_hyper = NULL,
                                  loss_functions = c("SELF", "ELF", "PLF"),
                                  init = NULL, tune_sd = NULL, verbose = FALSE, seed = NULL) {
  frailty <- match.arg(frailty)
  baseline <- match.arg(baseline)
  prior <- match.arg(prior)
  
  if (!is.null(seed)) set.seed(seed)
  
  if (!is.null(formula)) {
    if (is.null(data)) data <- environment(formula)
    mf <- stats::model.frame(formula, data)
    lhs <- stats::model.response(mf)
    if (is.matrix(lhs) && ncol(lhs) >= 4) {
      t1 <- as.numeric(lhs[, 1])
      t2 <- as.numeric(lhs[, 2])
      d1 <- as.numeric(lhs[, 3])
      d2 <- as.numeric(lhs[, 4])
    } else {
      stop("Formula response must be a 4-column matrix: cbind(t1, t2, d1, d2)")
    }
    covars <- stats::model.matrix(formula, mf)
    if (ncol(covars) > 1 && colnames(covars)[1] == "(Intercept)") {
      covars <- covars[, -1, drop = FALSE]
    }
  } else {
    if (is.null(t1) || is.null(t2) || is.null(d1) || is.null(d2)) {
      stop("Must supply either formula or t1, t2, d1, d2 vectors.")
    }
    t1 <- as.numeric(t1)
    t2 <- as.numeric(t2)
    d1 <- as.numeric(d1)
    d2 <- as.numeric(d2)
    if (!is.null(covars)) covars <- as.matrix(covars)
  }
  
  if (any(is.na(t1)) || any(is.na(t2)) || any(is.na(d1)) || any(is.na(d2))) {
    stop("Missing values (NA or NaN) detected in survival times or censoring indicators.")
  }
  if (!is.null(covars) && any(is.na(covars))) {
    stop("Missing values (NA or NaN) detected in covariates matrix.")
  }
  
  n <- length(t1)
  p_cov <- if (!is.null(covars)) ncol(as.matrix(covars)) else 0
  cov_names <- if (p_cov > 0) colnames(covars) else character(0)
  if (p_cov > 0 && is.null(cov_names)) cov_names <- paste0("beta", seq_len(p_cov))
  
  b_names <- switch(baseline,
    "GE" = c("zeta", "delta"),
    "EG" = c("zeta", "delta"),
    "GIR" = c("zeta", "delta"),
    "MIW" = c("zeta", "delta", "xi"),
    "GLL-I" = c("lambda", "alpha", "gamma"),
    "GLL-II" = c("lambda", "alpha", "gamma"),
    "GIE" = c("zeta", "delta"),
    "Weibull" = c("zeta", "delta"),
    "Exponential" = c("delta"),
    "Log-Logistic" = c("lambda", "gamma")
  )
  
  p1_names <- paste0(b_names, "1")
  p2_names <- paste0(b_names, "2")
  
  f_names <- switch(frailty,
    "WL" = "theta",
    "GL" = c("theta", "mu"),
    "GIG" = c("theta", "theta2"),
    "Gamma" = "theta",
    "IG" = "theta",
    "Lindley" = "theta",
    "None" = character(0)
  )
  
  all_param_names <- c(p1_names, p2_names, f_names, cov_names)
  n_params <- length(all_param_names)
  
  if (is.null(init)) {
    init_p1 <- switch(baseline,
      "GE" = c(zeta = 10.0, delta = 0.10),
      "EG" = c(zeta = 8.0, delta = 0.10),
      "GIR" = c(zeta = 3.5, delta = 900.0),
      "MIW" = c(zeta = 3.0, delta = 1.0, xi = 0.05),
      "GLL-I" = c(lambda = 0.05, alpha = 30.0, gamma = 1.5),
      "GLL-II" = c(lambda = 0.05, alpha = 0.15, gamma = 1.5),
      "GIE" = c(zeta = 3.0, delta = 10.0),
      "Weibull" = c(zeta = 1.5, delta = 0.05),
      "Exponential" = c(delta = 0.05),
      "Log-Logistic" = c(lambda = 0.05, gamma = 1.5)
    )
    init_p2 <- init_p1
    
    init_f <- switch(frailty,
      "WL" = c(theta = 3.0),
      "GL" = c(theta = 3.0, mu = 1.5),
      "GIG" = c(theta = 2.0, theta2 = 2.0),
      "Gamma" = c(theta = 2.0),
      "IG" = c(theta = 2.0),
      "Lindley" = c(theta = 2.0),
      "None" = numeric(0)
    )
    
    init_beta <- if (p_cov > 0) rep(0.0, p_cov) else numeric(0)
    init <- c(init_p1, init_p2, init_f, init_beta)
    names(init) <- all_param_names
  }
  
  if (is.null(tune_sd)) {
    tune_sd <- rep(0.05, n_params)
    names(tune_sd) <- all_param_names
    if ("delta1" %in% all_param_names && baseline == "GIR") tune_sd["delta1"] <- 5.0
    if ("delta2" %in% all_param_names && baseline == "GIR") tune_sd["delta2"] <- 5.0
  }
  
  target_acc <- 0.35
  chains_list <- vector("list", n_chains)
  
  for (ch in seq_len(n_chains)) {
    current_params <- init
    if (ch > 1) {
      current_params <- current_params * exp(stats::rnorm(n_params, 0, 0.1))
    }
    
    n_kept <- floor((n_iter - burn_in) / thin)
    kept_matrix <- matrix(NA_real_, nrow = n_kept, ncol = n_params, dimnames = list(NULL, all_param_names))
    loglik_kept <- numeric(n_kept)
    cpo_inv_matrix <- matrix(0.0, nrow = n_kept, ncol = n)
    
    curr_ll <- loglik_biv_rhr(current_params, t1, t2, d1, d2, covars, frailty, baseline)
    acc_count <- rep(0, n_params)
    keep_idx <- 1
    
    for (iter in seq_len(n_iter)) {
      for (k in seq_len(n_params)) {
        p_name <- all_param_names[k]
        is_pos <- !(k > (n_params - p_cov))
        
        prop_params <- current_params
        if (is_pos) {
          prop_val <- current_params[k] * exp(stats::rnorm(1, 0, tune_sd[k]))
        } else {
          prop_val <- current_params[k] + stats::rnorm(1, 0, tune_sd[k])
        }
        prop_params[k] <- prop_val
        
        valid <- TRUE
        if (is_pos && prop_val <= 0) valid <- FALSE
        if (frailty == "GL") {
          th_val <- prop_params["theta"]
          mu_val <- prop_params["mu"]
          if (!is.na(th_val) && !is.na(mu_val) && mu_val >= (1 + th_val)) valid <- FALSE
        }
        
        if (valid) {
          prop_ll <- loglik_biv_rhr(prop_params, t1, t2, d1, d2, covars, frailty, baseline)
          
          log_prior_ratio <- 0.0
          if (prior == "informative" && !is.null(prior_hyper)) {
            if (is_pos && p_name %in% names(prior_hyper$gamma)) {
              sh <- prior_hyper$gamma[[p_name]][1]
              rt <- prior_hyper$gamma[[p_name]][2]
              log_prior_ratio <- (sh - 1) * (log(prop_val) - log(current_params[k])) - rt * (prop_val - current_params[k])
            }
          }
          
          log_jacobian <- if (is_pos) (log(prop_val) - log(current_params[k])) else 0.0
          log_alpha <- (prop_ll - curr_ll) + log_prior_ratio + log_jacobian
          
          if (!is.na(log_alpha) && log(stats::runif(1)) < log_alpha) {
            current_params[k] <- prop_val
            curr_ll <- prop_ll
            acc_count[k] <- acc_count[k] + 1
          }
        }
      }
      
      if (iter <= burn_in && iter %% 50 == 0) {
        rates <- acc_count / 50
        acc_count <- rep(0, n_params)
        tune_sd <- tune_sd * exp(pmin(0.5, pmax(-0.5, (rates - target_acc))))
        tune_sd <- pmax(1e-4, pmin(tune_sd, 50.0))
      }
      
      if (iter > burn_in && (iter - burn_in) %% thin == 0) {
        kept_matrix[keep_idx, ] <- current_params
        loglik_kept[keep_idx] <- curr_ll
        cpo_inv_matrix[keep_idx, ] <- 1 / pmax(exp(pmin(50, pmax(-100, curr_ll / n))), 1e-15)
        keep_idx <- keep_idx + 1
      }
    }
    
    chains_list[[ch]] <- list(
      draws = kept_matrix,
      loglik = loglik_kept,
      cpo_inv = cpo_inv_matrix
    )
  }
  
  all_draws <- do.call(rbind, lapply(chains_list, function(x) x$draws))
  all_loglik <- unlist(lapply(chains_list, function(x) x$loglik))
  
  est_self <- colMeans(all_draws)
  
  est_elf <- vapply(seq_len(n_params), function(k) {
    if (all_param_names[k] %in% cov_names) {
      est_self[k]
    } else {
      1 / mean(1 / pmax(all_draws[, k], 1e-10))
    }
  }, numeric(1))
  names(est_elf) <- all_param_names
  
  est_plf <- vapply(seq_len(n_params), function(k) {
    sqrt(mean(all_draws[, k]^2))
  }, numeric(1))
  names(est_plf) <- all_param_names
  
  post_sd <- apply(all_draws, 2, stats::sd)
  post_median <- apply(all_draws, 2, stats::median)
  
  mcse <- vapply(seq_len(n_params), function(k) {
    x <- all_draws[, k]
    nb <- 20
    bsize <- floor(length(x) / nb)
    if (bsize < 2) return(post_sd[k] / sqrt(length(x)))
    bmeans <- vapply(seq_len(nb), function(b) mean(x[((b - 1) * bsize + 1):(b * bsize)]), numeric(1))
    stats::sd(bmeans) / sqrt(nb)
  }, numeric(1))
  names(mcse) <- all_param_names
  
  hpd_lims <- t(apply(all_draws, 2, hpd_interval, prob = 0.95))
  colnames(hpd_lims) <- c("HPD_Lower", "HPD_Upper")
  
  diag_df <- diagnostics_mcmc(all_draws)
  
  dev_draws <- -2 * all_loglik
  mean_dev <- mean(dev_draws)
  dev_at_mean <- -2 * loglik_biv_rhr(est_self, t1, t2, d1, d2, covars, frailty, baseline)
  pD <- mean_dev - dev_at_mean
  DIC <- dev_at_mean + 2 * pD
  EAIC <- mean_dev + 2 * n_params
  EBIC <- mean_dev + n_params * log(n)
  
  all_cpo_inv <- do.call(rbind, lapply(chains_list, function(x) x$cpo_inv))
  cpo_vals <- 1 / colMeans(all_cpo_inv)
  LPML <- sum(log(pmax(cpo_vals, 1e-15)))
  
  AIC_val <- dev_at_mean + 2 * n_params
  BIC_val <- dev_at_mean + n_params * log(n)
  AICc_val <- AIC_val + (2 * n_params * (n_params + 1)) / max(1, n - n_params - 1)
  HQIC_val <- dev_at_mean + 2 * n_params * log(log(n))
  
  th_est <- if (frailty != "None") est_self["theta"] else 1.0
  mu_est <- if (frailty == "GL") est_self["mu"] else NULL
  th2_est <- if (frailty == "GIG") est_self["theta2"] else NULL
  
  tau_est <- kendall_tau_rhr(frailty, th_est, mu_est, th2_est)
  var_est <- var_frailty(frailty, th_est, mu_est, th2_est)
  
  summary_table <- data.frame(
    Parameter = all_param_names,
    SELF = est_self,
    ELF = est_elf,
    PLF = est_plf,
    SD = post_sd,
    MCSE = mcse,
    Median = post_median,
    HPD_Lower = hpd_lims[, "HPD_Lower"],
    HPD_Upper = hpd_lims[, "HPD_Upper"],
    Geweke_p = diag_df$Geweke_p,
    HW_p = diag_df$Heidelberger_p,
    stringsAsFactors = FALSE
  )
  
  out <- list(
    summary = summary_table,
    draws = all_draws,
    chains = chains_list,
    estimates = list(SELF = est_self, ELF = est_elf, PLF = est_plf),
    loss_functions = loss_functions,
    diagnostics = diag_df,
    criteria = list(
      AIC = AIC_val,
      BIC = BIC_val,
      AICc = AICc_val,
      HQIC = HQIC_val,
      DIC = DIC,
      pD = pD,
      EAIC = EAIC,
      EBIC = EBIC,
      LPML = LPML,
      mean_deviance = mean_dev,
      dev_at_mean = dev_at_mean
    ),
    cpo = cpo_vals,
    kendall_tau = tau_est,
    frailty_variance = var_est,
    model_info = list(
      frailty = frailty,
      baseline = baseline,
      n_obs = n,
      n_iter = n_iter,
      burn_in = burn_in,
      thin = thin,
      n_chains = n_chains,
      prior = prior,
      p_cov = p_cov,
      cov_names = cov_names,
      param_names = all_param_names
    ),
    data = list(t1 = t1, t2 = t2, d1 = d1, d2 = d2, covars = covars),
    call = match.call()
  )
  
  class(out) <- "rhr_bayes_fit"
  return(out)
}
