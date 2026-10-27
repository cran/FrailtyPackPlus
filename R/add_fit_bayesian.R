#' Bayesian Estimation for Additive Shared Frailty Models via MCMC
#'
#' Fits additive shared frailty regression models for bivariate survival data
#' using an adaptive hybrid Metropolis-Hastings Markov Chain Monte Carlo (MCMC)
#' sampler with multiple chains, informative or flat priors, and point estimation
#' under Squared Error, Modified Quadratic Squared Error, and Precautionary loss functions.
#'
#' @param formula Model formula for covariates (e.g. \code{~ age + sex + gn + an + pkd}).
#'   If \code{NULL}, an intercept-only or single-covariate model is fitted from the data.
#' @param data A data frame containing survival times, censoring indicators, and covariates.
#' @param t1_col Character string naming the first survival time column in \code{data}.
#' @param delta1_col Character string naming the first censoring indicator column in \code{data}.
#' @param t2_col Character string naming the second survival time column in \code{data}.
#' @param delta2_col Character string naming the second censoring indicator column in \code{data}.
#' @param frailty_type Character string specifying the frailty distribution:
#'   \code{"genlindley"} (Generalized Lindley), \code{"wlindley"} (Weighted Lindley),
#'   \code{"gamma"} (Gamma), or \code{"invgaussian"} (Inverse Gaussian).
#' @param baseline_dist Character string specifying the baseline distribution:
#'   \code{"gw"} (Generalized Weibull) or \code{"gll2"} (Generalized Log-Logistic Type-II).
#' @param n_iter Total number of MCMC iterations per chain (default 10000).
#' @param burn_in Number of initial burn-in iterations to discard (default 2000).
#' @param thin Thinning interval / autocorrelation lag (default 5).
#' @param n_chains Number of parallel/independent Markov chains to run (default 2).
#' @param prior Character string specifying prior distribution: \code{"flat"} (non-informative)
#'   or \code{"informative"}.
#' @param prior_params Optional list of hyperparameters for informative priors.
#' @param init Optional list or vector of initial parameter values.
#' @param seed Optional integer seed for reproducibility.
#' @param verbose Logical; if \code{TRUE}, prints MCMC chain progress.
#'
#' @return An object of class \code{c("add_frailty_bayes", "add_frailty")} containing:
#' \describe{
#'   \item{estimates}{Data frame of Bayes point estimates (SELF, MQSELF, PLF), risks, and credible intervals.}
#'   \item{chains}{List of thinned posterior sample matrices for each chain.}
#'   \item{posterior_samples}{Combined posterior sample matrix across all chains.}
#'   \item{deviance_samples}{Posterior deviance samples (-2 * logLik).}
#'   \item{diagnostics}{List containing Geweke, Gelman-Rubin, and Heidelberger-Welch convergence statistics.}
#'   \item{criteria}{Model selection criteria (AIC, BIC, AICc, HQIC, DIC, LogLik).}
#'   \item{kendall_tau}{Kendall's tau dependence measure evaluated at the Bayes estimate.}
#'   \item{frailty_variance}{Theoretical frailty variance at the Bayes estimate.}
#'   \item{model_info}{List containing model settings, formula, call, and data dimensions.}
#' }
#'
#' @export
#' @references
#' Pandey, A., Hanagal, D. D., and Tyagi, S. (2022). Generalised Lindley shared
#' additive frailty regression model for bivariate survival data.
#' \emph{Statistics in Transition new series}, 23(4), 161-176.
#' \doi{10.2478/stattrans-2022-0048}
#'
#' Tyagi, S., Pandey, A., Hanagal, D. D., and Chesneau, C. (2024). A study on
#' comparisons of additive regression frailty models to counter heterogeneity:
#' Bayesian strategies and case study. \emph{Communications in Statistics -
#' Simulation and Computation}. \doi{10.1080/03610918.2024.2391874}
#'
#' @examples
#' \donttest{
#' data(kidney_bivariate)
#' sub_data <- kidney_bivariate[1:20, ]
#' fit_bayes <- fit_add_frailty_bayes(
#'   formula = ~ age,
#'   data = sub_data,
#'   frailty_type = "wlindley", baseline_dist = "gw",
#'   n_iter = 200, burn_in = 50, thin = 2, n_chains = 2,
#'   seed = 123
#' )
#' print(fit_bayes)
#' }
fit_add_frailty_bayes <- function(formula = NULL,
                                  data,
                                  t1_col = "t1",
                                  delta1_col = "delta1",
                                  t2_col = "t2",
                                  delta2_col = "delta2",
                                  frailty_type = c("genlindley", "wlindley", "gamma", "invgaussian"),
                                  baseline_dist = c("gw", "gll2"),
                                  n_iter = 10000,
                                  burn_in = 2000,
                                  thin = 5,
                                  n_chains = 2,
                                  prior = c("flat", "informative"),
                                  prior_params = NULL,
                                  init = NULL,
                                  seed = NULL,
                                  verbose = FALSE) {
  call <- match.call()
  if (!is.null(seed)) set.seed(seed)
  frailty_type <- match.arg(frailty_type)
  baseline_dist <- match.arg(baseline_dist)
  prior <- match.arg(prior)

  if (!is.data.frame(data)) {
    stop("'data' must be a data frame.")
  }

  t1 <- as.numeric(data[[t1_col]])
  delta1 <- as.numeric(data[[delta1_col]])
  t2 <- as.numeric(data[[t2_col]])
  delta2 <- as.numeric(data[[delta2_col]])

  if (is.null(formula)) {
    cov_cols <- setdiff(colnames(data), c(t1_col, delta1_col, t2_col, delta2_col, "id", "w"))
    if (length(cov_cols) > 0) {
      X <- as.matrix(data[, cov_cols, drop = FALSE])
    } else {
      X <- matrix(0, nrow = length(t1), ncol = 0)
    }
  } else {
    mf <- stats::model.frame(formula, data)
    X <- stats::model.matrix(formula, data)
    # Remove intercept if present, as additive hazard handles baseline
    if ("(Intercept)" %in% colnames(X)) {
      X <- X[, colnames(X) != "(Intercept)", drop = FALSE]
    }
  }

  n_cov <- ncol(X)
  cov_names <- if (n_cov > 0) colnames(X) else character(0)

  # Parameter names
  if (frailty_type == "genlindley") {
    par_names <- c("zeta1", "delta1", "xi1", "zeta2", "delta2", "xi2", "theta", "mu")
  } else if (frailty_type == "wlindley") {
    par_names <- c("zeta1", "delta1", "xi1", "zeta2", "delta2", "xi2", "s")
  } else {
    par_names <- c("zeta1", "delta1", "xi1", "zeta2", "delta2", "xi2", "theta")
  }
  if (n_cov > 0) {
    par_names <- c(par_names, cov_names)
  }
  n_pars <- length(par_names)
  n_pos_pars <- if (frailty_type == "genlindley") 8 else 7

  # Log-prior density function
  log_prior <- function(pars) {
    if (any(pars[1:n_pos_pars] <= 0)) return(-1e12)
    if (frailty_type == "genlindley") {
      theta_val <- pars[7]
      mu_val <- pars[8]
      if (mu_val >= 1 + theta_val) return(-1e12)
    }

    if (prior == "flat") {
      # Gamma(0.0001, 0.0001) for positive parameters
      lp_pos <- sum(stats::dgamma(pars[1:n_pos_pars], shape = 0.0001, rate = 0.0001, log = TRUE))
      # Normal(0, 1000^2) for regression coefficients
      lp_reg <- if (n_cov > 0) {
        sum(stats::dnorm(pars[(n_pos_pars + 1):n_pars], mean = 0, sd = 1000, log = TRUE))
      } else 0
      lp_pos + lp_reg
    } else {
      # Informative prior
      if (is.null(prior_params)) {
        # Default informative matching empirical moments
        lp_pos <- sum(stats::dgamma(pars[1:n_pos_pars], shape = 2, rate = 2 / pars[1:n_pos_pars], log = TRUE))
        lp_reg <- if (n_cov > 0) sum(stats::dnorm(pars[(n_pos_pars + 1):n_pars], mean = 0, sd = 10, log = TRUE)) else 0
        lp_pos + lp_reg
      } else {
        # User defined informative
        lp <- 0
        for (k in seq_len(n_pos_pars)) {
          shape_k <- if (!is.null(prior_params$shape[k])) prior_params$shape[k] else 2
          rate_k <- if (!is.null(prior_params$rate[k])) prior_params$rate[k] else 1
          lp <- lp + stats::dgamma(pars[k], shape = shape_k, rate = rate_k, log = TRUE)
        }
        if (n_cov > 0) {
          reg_mean <- if (!is.null(prior_params$beta_mean)) prior_params$beta_mean else rep(0, n_cov)
          reg_sd <- if (!is.null(prior_params$beta_sd)) prior_params$beta_sd else rep(10, n_cov)
          lp <- lp + sum(stats::dnorm(pars[(n_pos_pars + 1):n_pars], mean = reg_mean, sd = reg_sd, log = TRUE))
        }
        lp
      }
    }
  }

  log_posterior <- function(pars) {
    lp <- log_prior(pars)
    if (lp <= -1e10) return(-1e12)
    ll <- loglik_add_frailty(pars, t1, delta1, t2, delta2, X, frailty_type, baseline_dist)
    if (ll <= -1e10) return(-1e12)
    ll + lp
  }

  # Initial values generator
  get_init_pars <- function(chain_idx) {
    if (!is.null(init)) {
      if (is.list(init) && length(init) >= chain_idx) return(init[[chain_idx]])
      if (is.numeric(init) && length(init) == n_pars) return(init)
    }

    # Default robust initial points based on model
    if (baseline_dist == "gw") {
      base_init1 <- c(2.0 * (1 + 0.1 * (chain_idx - 1)), 0.05, 0.65)
      base_init2 <- c(2.7 * (1 + 0.1 * (chain_idx - 1)), 0.06, 0.67)
    } else {
      base_init1 <- c(3.8 * (1 + 0.1 * (chain_idx - 1)), 0.05, 1.2)
      base_init2 <- c(4.5 * (1 + 0.1 * (chain_idx - 1)), 0.06, 0.9)
    }

    if (frailty_type == "genlindley") {
      frail_init <- c(2.3 * (1 + 0.1 * (chain_idx - 1)), 1.4)
    } else if (frailty_type == "wlindley") {
      frail_init <- c(1.1 * (1 + 0.1 * (chain_idx - 1)))
    } else {
      frail_init <- c(1.0 * (1 + 0.1 * (chain_idx - 1)))
    }

    beta_init <- if (n_cov > 0) rep(0.01 * chain_idx, n_cov) else numeric(0)
    c(base_init1, base_init2, frail_init, beta_init)
  }

  chains_samples <- vector("list", n_chains)
  chains_deviances <- vector("list", n_chains)

  for (c_idx in seq_len(n_chains)) {
    if (verbose) cat(sprintf("Running MCMC chain %d/%d...\n", c_idx, n_chains))
    current_pars <- get_init_pars(c_idx)
    current_logpost <- log_posterior(current_pars)

    # If initial point invalid, adjust
    if (current_logpost <= -1e10) {
      current_pars[1:n_pos_pars] <- abs(current_pars[1:n_pos_pars]) + 0.1
      if (frailty_type == "genlindley" && current_pars[8] >= 1 + current_pars[7]) {
        current_pars[8] <- (1 + current_pars[7]) * 0.5
      }
      current_logpost <- log_posterior(current_pars)
    }

    prop_sd <- rep(0.05, n_pars)
    prop_sd[1:n_pos_pars] <- 0.04
    if (n_cov > 0) prop_sd[(n_pos_pars + 1):n_pars] <- 0.08

    # Storage for kept draws
    n_kept <- floor((n_iter - burn_in) / thin)
    kept_draws <- matrix(0, nrow = n_kept, ncol = n_pars)
    kept_dev <- numeric(n_kept)
    colnames(kept_draws) <- par_names

    kept_count <- 0
    accepted_count <- 0

    for (iter in seq_len(n_iter)) {
      # Component-wise / Block adaptive Metropolis-Hastings
      for (p_idx in seq_len(n_pars)) {
        prop_pars <- current_pars
        if (p_idx <= n_pos_pars) {
          # Log-scale random walk proposal
          curr_log <- log(current_pars[p_idx])
          prop_log <- curr_log + stats::rnorm(1, mean = 0, sd = prop_sd[p_idx])
          prop_pars[p_idx] <- exp(prop_log)
          hastings_ratio <- prop_pars[p_idx] / current_pars[p_idx]
        } else {
          # Standard normal proposal for beta
          prop_pars[p_idx] <- current_pars[p_idx] + stats::rnorm(1, mean = 0, sd = prop_sd[p_idx])
          hastings_ratio <- 1.0
        }

        # Check boundary condition for GL frailty
        if (frailty_type == "genlindley") {
          if (prop_pars[8] >= 1 + prop_pars[7]) next
        }

        prop_logpost <- log_posterior(prop_pars)
        if (is.nan(prop_logpost) || is.na(prop_logpost)) prop_logpost <- -1e12

        alpha_acc <- (prop_logpost - current_logpost) + log(hastings_ratio)
        if (!is.nan(alpha_acc) && log(stats::runif(1)) < alpha_acc) {
          current_pars <- prop_pars
          current_logpost <- prop_logpost
          if (iter <= burn_in) accepted_count <- accepted_count + 1
        }
      }

      # Adapt proposal scale during burn-in every 100 iterations
      if (iter <= burn_in && iter %% 100 == 0) {
        acc_rate <- accepted_count / (100 * n_pars)
        if (acc_rate < 0.20) {
          prop_sd <- prop_sd * 0.8
        } else if (acc_rate > 0.50) {
          prop_sd <- prop_sd * 1.2
        }
        accepted_count <- 0
      }

      # Record post-burn-in samples
      if (iter > burn_in && (iter - burn_in) %% thin == 0) {
        kept_count <- kept_count + 1
        kept_draws[kept_count, ] <- current_pars
        # Deviance: -2 * loglik
        ll_curr <- loglik_add_frailty(current_pars, t1, delta1, t2, delta2, X, frailty_type, baseline_dist)
        kept_dev[kept_count] <- -2 * ll_curr
      }
    }

    chains_samples[[c_idx]] <- kept_draws
    chains_deviances[[c_idx]] <- kept_dev
  }

  all_posterior <- do.call(rbind, chains_samples)
  all_deviances <- unlist(chains_deviances)

  # Loss function estimates (SELF, MQSELF, PLF)
  est_df <- calc_loss_estimates(all_posterior)

  # Convergence diagnostics
  geweke_df <- calc_geweke(all_posterior)
  gr_df <- calc_gelman_rubin(chains_samples)
  hw_df <- calc_heidelberger_welch(all_posterior)

  diag_combined <- data.frame(
    Parameter = par_names,
    Geweke_Z = geweke_df$Geweke_Z,
    Geweke_p = geweke_df$p_value,
    Gelman_Rubin_PSRF = gr_df$PSRF,
    Heidelberger_Welch_p = hw_df$p_value,
    Stationary = hw_df$Stationary,
    stringsAsFactors = FALSE
  )

  # Criteria
  post_mean_pars <- est_df$Est_SELF
  dev_at_mean <- -2 * loglik_add_frailty(post_mean_pars, t1, delta1, t2, delta2, X, frailty_type, baseline_dist)
  crit_list <- calc_model_criteria(
    loglik = -0.5 * dev_at_mean,
    k = n_pars,
    n = length(t1),
    deviance_samples = all_deviances,
    deviance_at_mean = dev_at_mean
  )

  # Dependence measures
  if (frailty_type == "genlindley") {
    kt <- kendall_tau("genlindley", theta = post_mean_pars[7], mu = post_mean_pars[8])
    fv <- frailty_variance("genlindley", theta = post_mean_pars[7], mu = post_mean_pars[8])
  } else if (frailty_type == "wlindley") {
    kt <- kendall_tau("wlindley", s_param = post_mean_pars[7])
    fv <- frailty_variance("wlindley", s_param = post_mean_pars[7])
  } else if (frailty_type == "gamma") {
    kt <- kendall_tau("gamma", theta = post_mean_pars[7])
    fv <- frailty_variance("gamma", theta = post_mean_pars[7])
  } else {
    kt <- kendall_tau("invgaussian", theta = post_mean_pars[7])
    fv <- frailty_variance("invgaussian", theta = post_mean_pars[7])
  }

  out <- list(
    estimates = est_df,
    coefficients = stats::setNames(post_mean_pars, par_names),
    chains = chains_samples,
    posterior_samples = all_posterior,
    deviance_samples = all_deviances,
    diagnostics = diag_combined,
    criteria = crit_list,
    kendall_tau = kt,
    frailty_variance = fv,
    model_info = list(
      frailty_type = frailty_type,
      baseline_dist = baseline_dist,
      prior = prior,
      n_obs = length(t1),
      n_pars = n_pars,
      par_names = par_names,
      formula = formula,
      call = call
    )
  )

  class(out) <- c("add_frailty_bayes", "add_frailty")
  out
}
