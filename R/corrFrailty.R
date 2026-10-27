#' Control Parameters for corrFrailty Estimation
#'
#' @description
#' Auxiliary function for setting MCMC and MLE optimization control parameters
#' in \code{\link{corrFrailty}}.
#'
#' @param n.iter Integer, total number of MCMC iterations per chain. Default is 5000.
#' @param n.burnin Integer, number of burn-in iterations to discard. Default is 1000.
#' @param n.thin Integer, thinning interval. Default is 5.
#' @param n.chains Integer, number of parallel/independent MCMC chains. Default is 2.
#' @param init Optional list of initial value vectors (one per chain).
#' @param tune Optional numeric vector of initial proposal standard deviations.
#' @param adapt.iter Integer, number of iterations during burn-in to adapt proposal standard deviations. Default is 500.
#' @param priors Optional list specifying prior hyperparameters.
#' @param seed Optional integer random seed for reproducibility.
#' @param optim.method Character, optimization algorithm for MLE. Default is \code{"BFGS"}.
#' @param optim.control List of control parameters passed to \code{\link[stats]{optim}}.
#'
#' @return A list with specified control parameters.
#'
#' @author
#' Shikhar Tyagi, Arvind Pandey, Vrijesh Tripathi
#'
#' @export
#' @examples
#' ctrl <- corrFrailty.control(n.iter = 1000, n.burnin = 200, n.chains = 2)
corrFrailty.control <- function(n.iter = 5000,
                                n.burnin = 1000,
                                n.thin = 5,
                                n.chains = 2,
                                init = NULL,
                                tune = NULL,
                                adapt.iter = 500,
                                priors = list(),
                                seed = NULL,
                                optim.method = "BFGS",
                                optim.control = list()) {
  list(
    n.iter = as.integer(n.iter),
    n.burnin = as.integer(n.burnin),
    n.thin = as.integer(n.thin),
    n.chains = as.integer(n.chains),
    init = init,
    tune = tune,
    adapt.iter = as.integer(adapt.iter),
    priors = priors,
    seed = seed,
    optim.method = optim.method,
    optim.control = optim.control
  )
}

#' Fit Correlated Frailty Models for Bivariate Survival Data
#'
#' @description
#' Fits parametric correlated frailty models to bivariate right-censored survival data
#' using Bayesian Markov Chain Monte Carlo (MCMC) simulation or Maximum Likelihood Estimation (MLE).
#' Supports Correlated Gamma, Inverse Gaussian, Positive Stable, and Compound Poisson frailty
#' distributions paired with Generalized Weibull, Generalized Log-Logistic, Weibull, Exponential,
#' Pareto, Linear Failure Rate, and Burr (Type XII) baseline hazards.
#'
#' @param formula A formula object. For wide format data, use \code{Surv(t1, d1, t2, d2) ~ x1 + x2}
#'   or \code{cbind(t1, d1, t2, d2) ~ x1 + x2}. For long format data with paired clusters,
#'   use \code{Surv(time, status) ~ x1 + x2} with argument \code{id}.
#' @param data A data.frame containing the survival times, censoring indicators, and covariates.
#' @param id Optional variable identifying paired observations in long format data (must have exactly 2 observations per cluster).
#' @param frailty Character specifying the frailty distribution. One of \code{"gamma"} (default),
#'   \code{"invgaussian"}, \code{"posstable"}, \code{"comppoisson"}, or \code{"none"}.
#' @param baseline Character specifying the baseline hazard distribution for the first survival time T1.
#'   One of \code{"genweibull"}, \code{"genloglogistic"}, \code{"weibull"} (default),
#'   \code{"exponential"}, \code{"pareto"}, \code{"lfr"}, or \code{"burr"}.
#' @param baseline2 Character specifying the baseline hazard distribution for the second survival time T2.
#'   Defaults to the same distribution as \code{baseline}.
#' @param method Character specifying the estimation method: \code{"mcmc"} (Bayesian MCMC, default)
#'   or \code{"mle"} (Maximum Likelihood Estimation).
#' @param control List of control parameters created by \code{\link{corrFrailty.control}}.
#' @param ... Additional arguments passed to control parameters.
#'
#' @return An object of class \code{"corrFrailty"} containing:
#' \item{coefficients}{Point estimates (posterior means for MCMC, MLEs for MLE).}
#' \item{summary.table}{Detailed table of estimates, standard errors, credible/confidence intervals, and diagnostics.}
#' \item{method}{Estimation method used (\code{"mcmc"} or \code{"mle"}).}
#' \item{frailty}{Frailty distribution.}
#' \item{baseline1}{Baseline hazard distribution for T1.}
#' \item{baseline2}{Baseline hazard distribution for T2.}
#' \item{logLik}{Log-likelihood evaluated at parameter estimates.}
#' \item{deviance}{Deviance (-2 * logLik).}
#' \item{AIC}{Akaike Information Criterion.}
#' \item{BIC}{Bayesian Information Criterion.}
#' \item{DIC}{Deviance Information Criterion (for MCMC).}
#' \item{pD}{Effective number of parameters (for MCMC).}
#' \item{LPML}{Log Pseudo-Marginal Likelihood (for MCMC).}
#' \item{chains}{List of MCMC chains (for MCMC).}
#' \item{all.samples}{Matrix of pooled posterior samples (for MCMC).}
#' \item{call}{The matched function call.}
#' \item{data}{The processed data matrix used for model fitting.}
#'
#' @references
#' Hanagal, D. D., Pandey, A., and Ganguly, A. (2017). Correlated gamma frailty
#' models for bivariate survival data. \emph{Communications in Statistics - Simulation
#' and Computation}, 46(5), 3627-3644. <doi:10.1080/03610918.2015.1085559>
#'
#' Hanagal, D. D. and Pandey, A. (2018). Correlated inverse Gaussian frailty
#' models for bivariate survival data. \emph{Communications in Statistics - Theory
#' and Methods}, 49(14), 3455-3475. <doi:10.1080/03610926.2018.1549256>
#'
#' Hanagal, D. D. (2020). Correlated positive stable frailty models.
#' \emph{Communications in Statistics - Theory and Methods}, 51(8), 2419-2436.
#' <doi:10.1080/03610926.2020.1736305>
#'
#' Hanagal, D. D. (2020). Analysis of Kidney infection Data Using Correlated
#' Inverse Gaussian Frailty Model. \emph{Statistics and Applications}, 18(1), 1-19.
#'
#' @author
#' Shikhar Tyagi, Arvind Pandey, Vrijesh Tripathi
#'
#' @export
#' @examples
#' # Simulate small bivariate survival dataset
#' set.seed(123)
#' dat <- sim_corrFrailty(n = 30, frailty = "gamma",
#'                        frailty.params = list(sigma = 0.5, rho = 0.4),
#'                        baseline = "exponential",
#'                        baseline.params = list(lambda = 0.05),
#'                        beta = c(0.3), seed = 123)
#'
#' # Fit Correlated Gamma Frailty Model using MCMC
#' fit <- corrFrailty(Surv(t1, d1, t2, d2) ~ x1, data = dat,
#'                    frailty = "gamma", baseline = "exponential",
#'                    method = "mcmc",
#'                    control = corrFrailty.control(n.iter = 100, n.burnin = 20,
#'                                                  n.thin = 1, n.chains = 1))
#' print(fit)
#' summary(fit)
corrFrailty <- function(formula,
                        data,
                        id = NULL,
                        frailty = c("gamma", "invgaussian", "posstable", "comppoisson", "none"),
                        baseline = c("genweibull", "genloglogistic", "weibull",
                                     "exponential", "pareto", "lfr", "burr"),
                        baseline2 = baseline,
                        method = c("mcmc", "mle"),
                        control = corrFrailty.control(),
                        ...) {
  cl <- match.call()
  frailty <- match.arg(frailty)
  baseline <- match.arg(baseline)
  baseline2 <- match.arg(baseline2)
  method <- match.arg(method)
  
  # Merge additional control arguments if passed in ...
  dots <- list(...)
  if (length(dots) > 0) {
    for (nm in names(dots)) {
      if (nm %in% names(control)) {
        control[[nm]] <- dots[[nm]]
      }
    }
  }
  
  if (missing(data) || is.null(data)) {
    stop("Argument 'data' must be supplied as a data.frame")
  }
  
  # Extract response and covariates
  lhs <- formula[[2]]
  rhs_terms <- stats::delete.response(stats::terms(formula))
  
  if (!is.null(substitute(id)) && as.character(substitute(id)) %in% names(data)) {
    # Long format paired data
    id_name <- as.character(substitute(id))
    id_vec <- data[[id_name]]
    unique_ids <- unique(id_vec)
    n_pairs <- length(unique_ids)
    
    # Evaluate 2-element response
    t_raw <- eval(lhs[[2]], data, parent.frame())
    d_raw <- as.integer(eval(lhs[[3]], data, parent.frame()))
    
    t1 <- numeric(n_pairs); d1 <- integer(n_pairs)
    t2 <- numeric(n_pairs); d2 <- integer(n_pairs)
    
    mf_rhs <- stats::model.frame(rhs_terms, data = data)
    X_raw <- stats::model.matrix(rhs_terms, data = mf_rhs)
    if ("(Intercept)" %in% colnames(X_raw)) {
      X_raw <- X_raw[, colnames(X_raw) != "(Intercept)", drop = FALSE]
    }
    
    X_mat <- matrix(0, nrow = n_pairs, ncol = ncol(X_raw))
    colnames(X_mat) <- colnames(X_raw)
    
    for (i in seq_along(unique_ids)) {
      cur_idx <- which(id_vec == unique_ids[i])
      if (length(cur_idx) != 2) {
        stop(paste("Each cluster in long format must have exactly 2 observations. Cluster",
                   unique_ids[i], "has", length(cur_idx)))
      }
      t1[i] <- t_raw[cur_idx[1]]; d1[i] <- d_raw[cur_idx[1]]
      t2[i] <- t_raw[cur_idx[2]]; d2[i] <- d_raw[cur_idx[2]]
      if (ncol(X_raw) > 0) {
        X_mat[i, ] <- colMeans(X_raw[cur_idx, , drop = FALSE])
      }
    }
    X <- X_mat
  } else if (is.call(lhs) && length(lhs) == 5) {
    # Wide format 4-element response: Surv(t1, d1, t2, d2) or cbind(t1, d1, t2, d2)
    t1 <- as.numeric(eval(lhs[[2]], data, parent.frame()))
    d1 <- as.integer(eval(lhs[[3]], data, parent.frame()))
    t2 <- as.numeric(eval(lhs[[4]], data, parent.frame()))
    d2 <- as.integer(eval(lhs[[5]], data, parent.frame()))
    
    mf_rhs <- stats::model.frame(rhs_terms, data = data)
    X <- stats::model.matrix(rhs_terms, data = mf_rhs)
    if ("(Intercept)" %in% colnames(X)) {
      X <- X[, colnames(X) != "(Intercept)", drop = FALSE]
    }
  } else {
    stop("Formula response must have 4 components: Surv(t1, d1, t2, d2) or specify 'id' for paired long format.")
  }
  
  # Data validation
  if (length(t1) != length(t2) || length(t1) != length(d1) || length(t1) != length(d2)) {
    stop("Lengths of t1, d1, t2, d2 must match.")
  }
  if (any(t1 <= 0) || any(t2 <= 0)) {
    stop("Survival times t1 and t2 must be strictly positive.")
  }
  if (any(!d1 %in% c(0, 1)) || any(!d2 %in% c(0, 1))) {
    stop("Censoring indicators d1 and d2 must be binary (0 or 1).")
  }
  
  if (!is.null(X) && ncol(X) == 0) X <- NULL
  
  # Fit Model
  if (method == "mcmc") {
    fit_res <- fit_corrFrailty_mcmc(t1, d1, t2, d2, X,
                                    baseline, baseline2, frailty,
                                    control)
  } else {
    fit_res <- fit_corrFrailty_mle(t1, d1, t2, d2, X,
                                   baseline, baseline2, frailty,
                                   control)
  }
  
  fit_res$call <- cl
  fit_res$formula <- formula
  fit_res$frailty <- frailty
  fit_res$baseline1 <- baseline
  fit_res$baseline2 <- baseline2
  fit_res$n.obs <- length(t1)
  fit_res$n.events1 <- sum(d1)
  fit_res$n.events2 <- sum(d2)
  fit_res$data <- list(t1 = t1, d1 = d1, t2 = t2, d2 = d2, X = X)
  
  class(fit_res) <- "corrFrailty"
  fit_res
}
