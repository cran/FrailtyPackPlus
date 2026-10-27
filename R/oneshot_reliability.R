#' Reliability Evaluation for k-out-of-M Systems
#'
#' @description
#' Evaluates the exact system reliability function \eqn{R_k(t)} of \eqn{k}-out-of-\eqn{M}
#' structured one-shot devices under shared gamma frailty at specified operating times
#' and stress levels.
#'
#' @param object A fitted model object of class \code{"oneshot_kofm"} or \code{"oneshot_frailty"}.
#' @param ... Additional arguments passed to methods.
#'
#' @return A data frame containing mission times, stress levels, system thresholds \eqn{k},
#'   and evaluated system reliabilities.
#' @export
reliability_k_out_of_m <- function(object, ...) {
  UseMethod("reliability_k_out_of_m")
}

#' Mean Time to Failure (MTTF) for k-out-of-M Systems
#'
#' @description
#' Computes the Mean Time to Failure (MTTF) of \eqn{k}-out-of-\eqn{M} systems under shared gamma frailty.
#'
#' @param object A fitted model object of class \code{"oneshot_kofm"} or \code{"oneshot_frailty"}.
#' @param ... Additional arguments passed to methods.
#'
#' @return A data frame containing MTTF values across system configurations.
#' @export
mttf_k_out_of_m <- function(object, ...) {
  UseMethod("mttf_k_out_of_m")
}

#' Expected Component Lifetimes and Moments
#'
#' @description
#' Calculates the expected mean failure times and moments for individual components under gamma frailty.
#'
#' @param object A fitted model object of class \code{"oneshot_kofm"} or \code{"oneshot_frailty"}.
#' @param ... Additional arguments passed to methods.
#'
#' @return A data frame or list containing component mean lifetimes and standard deviations.
#' @export
component_moments <- function(object, ...) {
  UseMethod("component_moments")
}


#' Reliability Evaluation for k-Out-of-M Systems
#'
#' @description
#' Evaluates the exact system reliability function \eqn{R_k(t)} of \eqn{k}-out-of-\eqn{M}
#' structured one-shot devices under shared gamma frailty at specified operating times
#' and stress levels.
#'
#' @param object A fitted model object of class \code{"oneshot_kofm"}.
#' @param t A positive numeric vector of mission times \eqn{t}.
#' @param stress A numeric vector or matrix specifying stress conditions.
#' @param k An optional integer vector specifying system structures (\eqn{1 \le k \le M}).
#'   If \code{NULL} (default), evaluates all structures \eqn{k = 1, \dots, M}.
#' @param ... Additional arguments passed to other methods.
#'
#' @return A data frame containing:
#'   \describe{
#'     \item{k}{The system structure threshold (\eqn{k}-out-of-\eqn{M}).}
#'     \item{structure}{Character description (\code{"Parallel (k=1)"}, \code{"Series (k=M)"}, or \code{"k-out-of-M"}).}
#'     \item{time}{Mission time \eqn{t}.}
#'     \item{stress}{Operating stress level.}
#'     \item{reliability}{Model-estimated system reliability \eqn{R_k(t)}.}
#'   }
#'
#' @details
#' For a \eqn{k}-out-of-\eqn{M} system, at least \eqn{k} components must function for the
#' device to work. Under gamma frailty, the system reliability is given by:
#' \deqn{R_k(t) = \sum_{n=k}^M \left( \sum_{d=0}^{n-k} (-1)^d \binom{n}{d} \right) \sum_{X \in \mathcal{X}_n} g(X, t)}
#' where \eqn{\mathcal{X}_n} is the set of all subsets of \eqn{\Omega} of size \eqn{n}, and
#' \eqn{g(X, t) = (1 + \beta \sum_{m \in X} (\lambda_m t)^{\alpha_m})^{-1/\beta}}.
#'
#' Missing values are not permitted.
#'
#' @references
#' Ling, M.-H., Balakrishnan, N., Yu, C., and So, H. Y. (2021). Inference for One-Shot
#' Devices with Dependent k-Out-of-M Structured Components under Gamma Frailty.
#' \emph{Mathematics}, 9(23), 3032. \doi{10.3390/math9233032}.
#'
#' @examples
#' data(classH)
#' fit <- oneshot_kofm(
#'   data = classH,
#'   time = "time",
#'   stress = "stress",
#'   counts = c("n_none", "n_T", "n_G", "n_both"),
#'   baseline = "exponential"
#' )
#' rel <- reliability_k_out_of_m(fit, t = c(1000, 2000, 5000), stress = 356)
#' head(rel)
#' @export
reliability_k_out_of_m.oneshot_kofm <- function(object, t, stress, k = NULL, ...) {
  if (!inherits(object, "oneshot_kofm")) {
    stop("'object' must be a fitted object of class 'oneshot_kofm'.")
  }
  
  if (missing(t) || is.null(t) || any(is.na(t)) || any(t < 0)) {
    stop("'t' must be a non-negative numeric vector.")
  }
  
  M <- object$M
  if (is.null(k)) {
    k <- seq_len(M)
  } else {
    if (any(is.na(k)) || any(k < 1) || any(k > M)) {
      stop("'k' must contain integers between 1 and M (", M, ").")
    }
    k <- as.integer(unique(k))
  }
  
  if (is.vector(stress)) {
    stress_mat <- matrix(stress, ncol = 1)
  } else {
    stress_mat <- as.matrix(stress)
  }
  
  X_cov <- cbind(1, stress_mat)
  p_cov <- ncol(X_cov)
  
  par <- object$coefficients
  idx <- 1
  a_params <- matrix(par[idx:(idx + p_cov * M - 1)], nrow = M, ncol = p_cov, byrow = TRUE)
  idx <- idx + p_cov * M
  
  if (object$baseline == "weibull") {
    b_params <- matrix(par[idx:(idx + p_cov * M - 1)], nrow = M, ncol = p_cov, byrow = TRUE)
    idx <- idx + p_cov * M
  } else {
    b_params <- NULL
  }
  
  beta <- if (object$frailty) par["beta"] else 0
  
  results <- list()
  row_count <- 1
  
  for (s_row in seq_len(nrow(X_cov))) {
    x_i <- matrix(X_cov[s_row, ], nrow = 1)
    s_val <- stress_mat[s_row, 1]
    
    lambda_vec <- numeric(M)
    alpha_vec <- numeric(M)
    for (m in seq_len(M)) {
      lambda_vec[m] <- exp(sum(x_i * a_params[m, ]))
      alpha_vec[m] <- if (object$baseline == "weibull") exp(sum(x_i * b_params[m, ])) else 1
    }
    
    for (t_val in t) {
      g_cache <- list()
      
      get_g <- function(subset_idx) {
        if (length(subset_idx) == 0) return(1)
        key <- paste(sort(subset_idx), collapse = ",")
        if (!is.null(g_cache[[key]])) return(g_cache[[key]])
        
        cum_h <- sum((lambda_vec[subset_idx] * t_val)^(alpha_vec[subset_idx]))
        val <- if (beta > 1e-10) (1 + beta * cum_h)^(-1 / beta) else exp(-cum_h)
        g_cache[[key]] <<- val
        val
      }
      
      for (k_val in k) {
        rel_val <- 0
        for (n in k_val:M) {
          # Binomial sum: sum_{d=0}^{n-k} (-1)^d * choose(n, d)
          d_sum <- sum(vapply(0:(n - k_val), function(d) (-1)^d * choose(n, d), numeric(1)))
          
          # Subsets of size n
          combs <- utils::combn(M, n, simplify = FALSE)
          subset_g_sum <- 0
          for (X_set in combs) {
            subset_g_sum <- subset_g_sum + get_g(X_set)
          }
          rel_val <- rel_val + d_sum * subset_g_sum
        }
        
        rel_val <- max(0, min(1, rel_val))
        
        struct_label <- if (k_val == 1) {
          "Parallel (k=1)"
        } else if (k_val == M) {
          "Series (k=M)"
        } else {
          paste0(k_val, "-out-of-", M)
        }
        
        results[[row_count]] <- data.frame(
          k = k_val,
          structure = struct_label,
          time = t_val,
          stress = s_val,
          reliability = rel_val,
          stringsAsFactors = FALSE
        )
        row_count <- row_count + 1
      }
    }
  }
  
  do.call(rbind, results)
}

#' Mean Lifetime and Confidence Intervals for k-Out-of-M Systems
#'
#' @description
#' Computes the exact analytical mean lifetime \eqn{\mu_{0,k}} (Mean Time to Failure, MTTF)
#' of \eqn{k}-out-of-\eqn{M} structured one-shot devices under gamma frailty at normal
#' operating stress \eqn{s_0}, standard errors via the delta method, and both Asymptotic
#' Confidence Intervals (ACI) and Log-Transformed Confidence Intervals (TCI).
#'
#' @param object A fitted model object of class \code{"oneshot_kofm"}.
#' @param stress A numeric scalar or vector representing the normal operating condition \eqn{s_0}.
#' @param k An optional integer vector specifying system structures (\eqn{1 \le k \le M}).
#'   If \code{NULL} (default), evaluates all structures \eqn{k = 1, \dots, M}.
#' @param ci_level A numeric scalar between 0 and 1 specifying the confidence level (default 0.95).
#' @param upper_limit Upper numerical integration limit for Weibull baseline (default 1e6).
#' @param ... Additional arguments passed to other methods.
#'
#' @return A data frame containing:
#'   \describe{
#'     \item{k}{The system structure threshold (\eqn{k}-out-of-\eqn{M}).}
#'     \item{structure}{Character description (\code{"Parallel (k=1)"}, \code{"Series (k=M)"}, or \code{"k-out-of-M"}).}
#'     \item{stress}{Operating stress level \eqn{s_0}.}
#'     \item{mttf}{Estimated mean lifetime \eqn{\mu_{0,k}}.}
#'     \item{se}{Standard error of mean lifetime \eqn{se(\mu)} from the delta method.}
#'     \item{aci_lower}{Lower bound of the Asymptotic Confidence Interval (ACI).}
#'     \item{aci_upper}{Upper bound of the Asymptotic Confidence Interval (ACI).}
#'     \item{tci_lower}{Lower bound of the Log-Transformed Confidence Interval (TCI).}
#'     \item{tci_upper}{Upper bound of the Log-Transformed Confidence Interval (TCI).}
#'   }
#'
#' @details
#' For Exponential baseline distributions under gamma frailty, the exact analytical
#' mean lifetime at stress level \eqn{s_0} is given by Equation (4) of Ling et al. (2021):
#' \deqn{\mu_{0,k} = \frac{1}{1-\beta} \sum_{n=k}^M \left( \sum_{d=0}^{n-k} (-1)^d \binom{n}{d} \right) \sum_{X \in \mathcal{X}_n} \left( \sum_{m \in X} \lambda_{0,m} \right)^{-1}}
#'
#' The standard error is computed via the delta method:
#' \deqn{se(\mu) = \sqrt{P_\mu V(\hat{\theta}) P_\mu^T}}
#' where \eqn{P_\mu = \nabla_\theta \mu_{0,k}} is the analytical gradient vector
#' (Section 5, Ling et al., 2021).
#'
#' The \eqn{100(1-\delta)\%} Log-Transformed Confidence Interval (TCI) is:
#' \deqn{\left( \hat{\mu}_{0,k} \exp\left( -\frac{z_{1-\delta/2} se(\mu)}{\hat{\mu}_{0,k}} \right), \; \hat{\mu}_{0,k} \exp\left( \frac{z_{1-\delta/2} se(\mu)}{\hat{\mu}_{0,k}} \right) \right)}
#' ensuring strictly positive lower bounds.
#'
#' Missing values are not permitted.
#'
#' @references
#' Ling, M.-H., Balakrishnan, N., Yu, C., and So, H. Y. (2021). Inference for One-Shot
#' Devices with Dependent k-Out-of-M Structured Components under Gamma Frailty.
#' \emph{Mathematics}, 9(23), 3032. \doi{10.3390/math9233032}.
#'
#' Balakrishnan, N. and Ling, M. H. (2014). Gamma lifetimes and one-shot device testing
#' analysis. \emph{Reliability Engineering & System Safety}, 126, 54--64.
#' \doi{10.1016/j.ress.2014.01.009}.
#'
#' @examples
#' data(classH)
#' fit <- oneshot_kofm(
#'   data = classH,
#'   time = "time",
#'   stress = "stress",
#'   counts = c("n_none", "n_T", "n_G", "n_both"),
#'   baseline = "exponential"
#' )
#' # Mean lifetimes at normal operating condition s0 = 356 F
#' mttf_k_out_of_m(fit, stress = 356)
#' @export
mttf_k_out_of_m.oneshot_kofm <- function(object, stress, k = NULL, ci_level = 0.95, upper_limit = 1e6, ...) {
  if (!inherits(object, "oneshot_kofm")) {
    stop("'object' must be a fitted object of class 'oneshot_kofm'.")
  }
  
  if (missing(stress) || is.null(stress) || any(is.na(stress))) {
    stop("'stress' must be a numeric vector.")
  }
  
  if (is.null(ci_level) || length(ci_level) != 1 || is.na(ci_level) || ci_level <= 0 || ci_level >= 1) {
    stop("'ci_level' must be a numeric scalar between 0 and 1.")
  }
  
  M <- object$M
  if (is.null(k)) {
    k <- seq_len(M)
  } else {
    if (any(is.na(k)) || any(k < 1) || any(k > M)) {
      stop("'k' must contain integers between 1 and M (", M, ").")
    }
    k <- as.integer(unique(k))
  }
  
  z_crit <- stats::qnorm(1 - (1 - ci_level) / 2)
  par <- object$coefficients
  vcov_mat <- object$vcov
  n_par <- length(par)
  
  s_vec <- as.numeric(stress)
  x0 <- c(1, s_vec)
  p_cov <- length(x0)
  
  results <- list()
  row_count <- 1
  
  for (s_val in s_vec) {
    x_i <- c(1, s_val)
    
    for (k_val in k) {
      # Function to compute mu_{0,k} given parameter vector theta
      calc_mu <- function(p_vec) {
        idx <- 1
        a_m <- matrix(p_vec[idx:(idx + p_cov * M - 1)], nrow = M, ncol = p_cov, byrow = TRUE)
        idx <- idx + p_cov * M
        
        if (object$baseline == "weibull") {
          b_m <- matrix(p_vec[idx:(idx + p_cov * M - 1)], nrow = M, ncol = p_cov, byrow = TRUE)
          idx <- idx + p_cov * M
        } else {
          b_m <- NULL
        }
        
        b_val <- if (object$frailty) p_vec["beta"] else 0
        if (b_val >= 1) b_val <- 0.999
        
        lam0 <- numeric(M)
        for (m in seq_len(M)) {
          lam0[m] <- exp(sum(x_i * a_m[m, ]))
        }
        
        if (object$baseline == "exponential") {
          # Analytical formula (Equation 4)
          total_mu <- 0
          for (n in k_val:M) {
            d_sum <- sum(vapply(0:(n - k_val), function(d) (-1)^d * choose(n, d), numeric(1)))
            combs <- utils::combn(M, n, simplify = FALSE)
            subset_sum <- 0
            for (X_set in combs) {
              subset_sum <- subset_sum + 1 / sum(lam0[X_set])
            }
            total_mu <- total_mu + d_sum * subset_sum
          }
          if (object$frailty && b_val > 0) {
            total_mu <- total_mu / (1 - b_val)
          }
          total_mu
        } else {
          # Numerical integration for Weibull baseline
          rel_fun <- function(t_grid) {
            vapply(t_grid, function(tt) {
              temp_obj <- object
              temp_obj$coefficients <- p_vec
              rel_df <- reliability_k_out_of_m(temp_obj, t = tt, stress = s_val, k = k_val)
              rel_df$reliability[1]
            }, numeric(1))
          }
          stats::integrate(rel_fun, lower = 0, upper = upper_limit, subdivisions = 500L)$value
        }
      }
      
      mu_est <- calc_mu(par)
      
      # Delta method gradient
      grad_mu <- tryCatch({
        numDeriv::grad(calc_mu, par)
      }, error = function(e) {
        rep(0, n_par)
      })
      
      var_mu <- as.numeric(t(grad_mu) %*% vcov_mat %*% grad_mu)
      se_mu <- sqrt(max(var_mu, 0))
      
      # ACI
      aci_l <- max(0, mu_est - z_crit * se_mu)
      aci_u <- mu_est + z_crit * se_mu
      
      # TCI (Log-Transformed)
      if (mu_est > 0 && se_mu > 0) {
        tci_factor <- exp(z_crit * se_mu / mu_est)
        tci_l <- mu_est / tci_factor
        tci_u <- mu_est * tci_factor
      } else {
        tci_l <- aci_l
        tci_u <- aci_u
      }
      
      struct_label <- if (k_val == 1) {
        "Parallel (k=1)"
      } else if (k_val == M) {
        "Series (k=M)"
      } else {
        paste0(k_val, "-out-of-", M)
      }
      
      results[[row_count]] <- data.frame(
        k = k_val,
        structure = struct_label,
        stress = s_val,
        mttf = mu_est,
        se = se_mu,
        aci_lower = aci_l,
        aci_upper = aci_u,
        tci_lower = tci_l,
        tci_upper = tci_u,
        stringsAsFactors = FALSE
      )
      row_count <- row_count + 1
    }
  }
  
  do.call(rbind, results)
}

#' Marginal Component Moments
#'
#' @description
#' Computes the marginal scale, mean, and variance for each failure mode under
#' Lomax / Exponential distributions with gamma frailty.
#'
#' @param object A fitted model object of class \code{"oneshot_kofm"}.
#' @param stress A numeric scalar or vector representing operating conditions.
#' @param ... Additional arguments passed to other methods.
#'
#' @return A data frame containing:
#'   \describe{
#'     \item{component}{Component name.}
#'     \item{stress}{Operating stress level.}
#'     \item{lambda}{Estimated scale parameter \eqn{\lambda_{0,m}}.}
#'     \item{mean}{Marginal mean lifetime \eqn{\mu_{T_m}}.}
#'     \item{variance}{Marginal lifetime variance \eqn{\sigma^2_{T_m}} (if \eqn{\beta < 0.5}).}
#'     \item{sd}{Marginal standard deviation.}
#'   }
#'
#' @details
#' Under gamma frailty, the marginal distribution of component \eqn{m} is Lomax
#' (Pareto Type II) with:
#' \deqn{\mu_{T_m} = \frac{1}{\lambda_m (1 - \beta)}, \quad \beta < 1}
#' \deqn{\sigma^2_{T_m} = \frac{1}{\lambda_m^2 (1 - \beta)^2 (1 - 2\beta)}, \quad \beta < 0.5}
#'
#' Missing values are not permitted.
#'
#' @references
#' Ling, M.-H., Balakrishnan, N., Yu, C., and So, H. Y. (2021). Inference for One-Shot
#' Devices with Dependent k-Out-of-M Structured Components under Gamma Frailty.
#' \emph{Mathematics}, 9(23), 3032. \doi{10.3390/math9233032}.
#'
#' @examples
#' data(classH)
#' fit <- oneshot_kofm(
#'   data = classH,
#'   time = "time",
#'   stress = "stress",
#'   counts = c("n_none", "n_T", "n_G", "n_both"),
#'   baseline = "exponential"
#' )
#' component_moments(fit, stress = 356)
#' @export
component_moments.oneshot_kofm <- function(object, stress, ...) {
  if (!inherits(object, "oneshot_kofm")) {
    stop("'object' must be a fitted object of class 'oneshot_kofm'.")
  }
  
  if (missing(stress) || is.null(stress) || any(is.na(stress))) {
    stop("'stress' must be a numeric vector.")
  }
  
  M <- object$M
  par <- object$coefficients
  beta <- if (object$frailty) par["beta"] else 0
  
  s_vec <- as.numeric(stress)
  x0 <- c(1, s_vec)
  p_cov <- length(x0)
  
  idx <- 1
  a_params <- matrix(par[idx:(idx + p_cov * M - 1)], nrow = M, ncol = p_cov, byrow = TRUE)
  
  results <- list()
  row_count <- 1
  
  for (s_val in s_vec) {
    x_i <- c(1, s_val)
    for (m in seq_len(M)) {
      lam <- exp(sum(x_i * a_params[m, ]))
      
      mean_val <- if (beta < 1) 1 / (lam * (1 - beta)) else Inf
      var_val <- if (beta < 0.5) 1 / (lam^2 * (1 - beta)^2 * (1 - 2 * beta)) else Inf
      sd_val <- sqrt(var_val)
      
      results[[row_count]] <- data.frame(
        component = object$components[m],
        stress = s_val,
        lambda = lam,
        mean = mean_val,
        variance = var_val,
        sd = sd_val,
        stringsAsFactors = FALSE
      )
      row_count <- row_count + 1
    }
  }
  
  do.call(rbind, results)
}


#' System Reliability for k-out-of-M Systems under Gamma Frailty
#'
#' @description
#' Computes the system reliability \eqn{R_k(t)} at specified mission times \eqn{t}
#' for \eqn{k}-out-of-\eqn{M} structured one-shot devices under shared gamma frailty
#' and Weibull or Exponential baseline lifetime distributions.
#'
#' @param object A fitted model object of class \code{"oneshot_frailty"}.
#' @param t A non-negative numeric vector of mission time points at which reliability is evaluated.
#' @param stress An optional numeric vector or matrix of stress covariate values.
#' @param k An integer specifying the system threshold (the system functions if at least \eqn{k} of its \eqn{M} components function).
#' @param ... Additional arguments passed to other methods.
#'
#' @return A data frame containing:
#'   \describe{
#'     \item{time}{Evaluation time point \eqn{t}.}
#'     \item{k}{System threshold \eqn{k}.}
#'     \item{reliability}{Computed system reliability \eqn{R_k(t)}.}
#'     \item{unreliability}{Failure probability \eqn{F_k(t) = 1 - R_k(t)}.}
#'   }
#'
#' @details
#' For a \eqn{k}-out-of-\eqn{M} system, the system operates if at least \eqn{k}
#' components survive. Using the inclusion-exclusion principle and the power set
#' \eqn{P(Omega)}, the unconditional system reliability is
#' \deqn{R_k(t) = \sum_{X: |X| \le M - k} P(X, t)}
#' where \eqn{|X|} is the number of malfunctioned components.
#'
#' Special cases include:
#' \itemize{
#'   \item \code{k = 1}: Parallel system (functions if at least 1 component functions).
#'   \item \code{k = M}: Series system (functions only if all \eqn{M} components function).
#' }
#'
#' Missing values are not permitted in \code{t} or \code{stress}.
#'
#' @references
#' Yu, C., So, H. Y., and Balakrishnan, N. (2025). \emph{Inference for one-shot devices
#' with Weibull component lifetimes under gamma frailty model}. In: Stochastic Modeling
#' and Statistical Methods, Elsevier Inc., Chapter 1, pp. 1--34.
#' \doi{10.1016/B978-0-44-331694-4.00006-2}.
#'
#' Ling, M. H., Balakrishnan, N., Yu, C., and So, H. Y. (2021). \emph{Inference for
#' one-shot devices with dependent k-out-of-M structured components under gamma frailty}.
#' Mathematics, 9(23), 3032. \doi{10.3390/math9233032}.
#'
#' @examples
#' data(classH)
#' start_par <- c(-3.997, -0.009, -14.412, 0.003, 0.004, -0.001, -5.493e-4, -0.005, 0.082)
#' fit_w <- oneshot_frailty(
#'   data = classH,
#'   time = "time",
#'   stress = "stress",
#'   counts = c("n_none", "n_T", "n_G", "n_both"),
#'   baseline = "weibull",
#'   frailty = TRUE,
#'   start = start_par
#' )
#' # Reliability of 1-out-of-2 (parallel) and 2-out-of-2 (series) systems at 400 deg F
#' reliability_k_out_of_m(fit_w, t = c(1000, 2000, 5000), stress = 400, k = 1)
#' reliability_k_out_of_m(fit_w, t = c(1000, 2000, 5000), stress = 400, k = 2)
#' @export
reliability_k_out_of_m.oneshot_frailty <- function(object, t, stress = NULL, k = NULL, ...) {
  if (!inherits(object, "oneshot_frailty")) {
    stop("'object' must be a fitted 'oneshot_frailty' model.")
  }
  if (missing(t) || is.null(t) || any(is.na(t)) || any(t < 0)) {
    stop("'t' must be a non-negative numeric vector.")
  }
  
  M <- object$M
  if (is.null(k)) {
    k_vals <- seq_len(M)
  } else {
    if (any(k < 1 | k > M | k != as.integer(k))) {
      stop(sprintf("'k' must be integer(s) between 1 and M (%d).", M))
    }
    k_vals <- as.integer(k)
  }
  
  p_cov <- ncol(object$X_cov)
  if (is.null(stress)) {
    if (p_cov == 1) {
      X_eval <- matrix(1, nrow = 1, ncol = 1)
    } else {
      X_eval <- matrix(colMeans(object$X_cov), nrow = 1)
    }
  } else {
    if (is.vector(stress)) {
      if (length(stress) == p_cov - 1) {
        X_eval <- matrix(c(1, stress), nrow = 1)
      } else if (length(stress) == p_cov && stress[1] == 1) {
        X_eval <- matrix(stress, nrow = 1)
      } else if (p_cov == 1) {
        X_eval <- matrix(1, nrow = 1, ncol = 1)
      } else {
        stop(sprintf("'stress' must have length %d (excluding intercept).", p_cov - 1))
      }
    } else if (is.matrix(stress) || is.data.frame(stress)) {
      X_eval <- as.matrix(stress)
      if (ncol(X_eval) == p_cov - 1) {
        X_eval <- cbind(1, X_eval)
      }
    }
  }
  
  powerset <- object$powerset %||% generate_powerset(M)
  subset_sizes <- vapply(powerset, length, integer(1))
  
  out_list <- list()
  row_id <- 1
  
  for (t_val in t) {
    prob_res <- compute_probabilities(object$coefficients, X_eval, t_val, M,
                                      object$baseline, object$frailty, powerset)
    p_vec <- prob_res$prob_matrix[1, ]
    
    for (k_curr in k_vals) {
      valid_subsets <- which(subset_sizes <= (M - k_curr))
      rel_val <- sum(p_vec[valid_subsets])
      rel_val <- min(max(rel_val, 0), 1)
      
      out_list[[row_id]] <- data.frame(
        time = t_val,
        k = k_curr,
        reliability = rel_val,
        unreliability = 1 - rel_val
      )
      row_id <- row_id + 1
    }
  }
  
  do.call(rbind, out_list)
}

#' Mean Time to Failure (MTTF) for k-out-of-M Systems
#'
#' @description
#' Calculates the Mean Time to Failure (MTTF) of \eqn{k}-out-of-\eqn{M} structured
#' one-shot devices by numerical integration of the system reliability function \eqn{R_k(t)}.
#'
#' @param object A fitted model object of class \code{"oneshot_frailty"}.
#' @param stress Optional numeric vector or matrix of stress covariate values.
#' @param ... Additional arguments passed to other methods.
#' @param k Optional integer specifying system threshold (defaults to all \eqn{k = 1, \dots, M}).
#' @param upper_limit Upper limit of integration in time units (default 1e5).
#' @param ... Additional arguments passed to other methods.
#'
#' @return A data frame containing:
#'   \describe{
#'     \item{k}{System threshold \eqn{k}.}
#'     \item{system_type}{Descriptive string (e.g. \code{"Parallel (k=1)"}, \code{"Series (k=M)"}, \code{"k-out-of-M"}).}
#'     \item{mttf}{Computed Mean Time to Failure.}
#'   }
#'
#' @details
#' The MTTF of a \eqn{k}-out-of-\eqn{M} system is given by
#' \deqn{MTTF_k = \int_0^{\infty} R_k(t)   dt}
#' evaluated using adaptive quadrature via \code{stats::integrate}.
#'
#' @references
#' Yu, C., So, H. Y., and Balakrishnan, N. (2025). \emph{Inference for one-shot devices
#' with Weibull component lifetimes under gamma frailty model}. In: Stochastic Modeling
#' and Statistical Methods, Elsevier Inc., Chapter 1, pp. 1--34.
#' \doi{10.1016/B978-0-44-331694-4.00006-2}.
#'
#' @examples
#' data(classH)
#' start_par <- c(-3.997, -0.009, -14.412, 0.003, 0.004, -0.001, -5.493e-4, -0.005, 0.082)
#' fit_w <- oneshot_frailty(
#'   data = classH,
#'   time = "time",
#'   stress = "stress",
#'   counts = c("n_none", "n_T", "n_G", "n_both"),
#'   baseline = "weibull",
#'   frailty = TRUE,
#'   start = start_par
#' )
#' mttf_k_out_of_m(fit_w, stress = 400, upper_limit = 50000)
#' @export
mttf_k_out_of_m.oneshot_frailty <- function(object, stress = NULL, k = NULL, upper_limit = 1e5, ...) {
  if (!inherits(object, "oneshot_frailty")) {
    stop("'object' must be a fitted 'oneshot_frailty' model.")
  }
  
  M <- object$M
  if (is.null(k)) {
    k_vals <- seq_len(M)
  } else {
    k_vals <- as.integer(k)
  }
  
  res_list <- list()
  for (i in seq_along(k_vals)) {
    k_curr <- k_vals[i]
    r_fun <- function(t_vec) {
      vapply(t_vec, function(tt) {
        reliability_k_out_of_m(object, t = tt, stress = stress, k = k_curr)$reliability
      }, numeric(1))
    }
    
    int_res <- stats::integrate(r_fun, lower = 0, upper = upper_limit, subdivisions = 500, rel.tol = 1e-4)
    
    sys_type <- if (k_curr == 1) {
      "Parallel (k=1)"
    } else if (k_curr == M) {
      sprintf("Series (k=%d)", M)
    } else {
      sprintf("%d-out-of-%d", k_curr, M)
    }
    
    res_list[[i]] <- data.frame(
      k = k_curr,
      system_type = sys_type,
      mttf = int_res$value
    )
  }
  
  do.call(rbind, res_list)
}

#' Marginal Moments for Components under Shared Gamma Frailty
#'
#' @description
#' Computes the marginal scale, shape, mean lifetime, variance, and standard deviation
#' for individual components under the Power Lomax marginal distribution.
#'
#' @param object A fitted model object of class \code{"oneshot_frailty"}.
#' @param stress Optional numeric vector or matrix of stress covariate values.
#' @param ... Additional arguments passed to other methods.
#'
#' @return A data frame containing:
#'   \describe{
#'     \item{component}{Component name.}
#'     \item{lambda}{Estimated scale parameter \eqn{\lambda}.}
#'     \item{alpha}{Estimated shape parameter \eqn{\alpha}.}
#'     \item{beta}{Frailty variance parameter \eqn{\beta}.}
#'     \item{mean}{Marginal mean lifetime (or \code{Inf} if \eqn{\beta \ge \alpha}).}
#'     \item{variance}{Marginal variance (or \code{Inf} if \eqn{\beta \ge \alpha/2}).}
#'     \item{sd}{Marginal standard deviation.}
#'   }
#'
#' @details
#' Under Weibull baseline lifetime and shared gamma frailty with variance \eqn{\beta},
#' the marginal unconditional lifetime distribution of component \eqn{m} is Power Lomax
#' with survival function
#' \deqn{S(t) = ( 1 + \beta (\lambda t)^\alpha )^{-1/\beta}}
#'
#' The \eqn{r}-th raw moment is given in closed form by:
#' \deqn{E[T^r] = \frac{\Gamma(1 + r/\alpha)   \Gamma(1/\beta - r/\alpha)}{\lambda^r   \beta^{r/\alpha}   \Gamma(1/\beta)}}
#' valid when \eqn{\beta < \alpha / r}.
#'
#' @references
#' Yu, C., So, H. Y., and Balakrishnan, N. (2025). \emph{Inference for one-shot devices
#' with Weibull component lifetimes under gamma frailty model}. In: Stochastic Modeling
#' and Statistical Methods, Elsevier Inc., Chapter 1, pp. 1--34.
#' \doi{10.1016/B978-0-44-331694-4.00006-2}.
#'
#' @examples
#' data(classH)
#' start_par <- c(-3.997, -0.009, -14.412, 0.003, 0.004, -0.001, -5.493e-4, -0.005, 0.082)
#' fit_w <- oneshot_frailty(
#'   data = classH,
#'   time = "time",
#'   stress = "stress",
#'   counts = c("n_none", "n_T", "n_G", "n_both"),
#'   baseline = "weibull",
#'   frailty = TRUE,
#'   start = start_par
#' )
#' component_moments(fit_w, stress = 400)
#' @export
component_moments.oneshot_frailty <- function(object, stress = NULL, ...) {
  if (!inherits(object, "oneshot_frailty")) {
    stop("'object' must be a fitted 'oneshot_frailty' model.")
  }
  
  M <- object$M
  p_cov <- ncol(object$X_cov)
  if (is.null(stress)) {
    if (p_cov == 1) {
      X_eval <- matrix(1, nrow = 1, ncol = 1)
    } else {
      X_eval <- matrix(colMeans(object$X_cov), nrow = 1)
    }
  } else {
    if (is.vector(stress)) {
      if (length(stress) == p_cov - 1) {
        X_eval <- matrix(c(1, stress), nrow = 1)
      } else if (length(stress) == p_cov && stress[1] == 1) {
        X_eval <- matrix(stress, nrow = 1)
      } else if (p_cov == 1) {
        X_eval <- matrix(1, nrow = 1, ncol = 1)
      } else {
        stop(sprintf("'stress' must have length %d (excluding intercept).", p_cov - 1))
      }
    } else {
      X_eval <- as.matrix(stress)
      if (ncol(X_eval) == p_cov - 1) {
        X_eval <- cbind(1, X_eval)
      }
    }
  }
  
  prob_res <- compute_probabilities(object$coefficients, X_eval, tau = 1, M = M,
                                    baseline = object$baseline, frailty = object$frailty)
  lambda_vec <- prob_res$lambda_mat[1, ]
  alpha_vec <- prob_res$alpha_mat[1, ]
  beta_val <- prob_res$beta
  
  mean_vec <- rep(NA_real_, M)
  var_vec <- rep(NA_real_, M)
  
  for (m in seq_len(M)) {
    lam <- lambda_vec[m]
    alp <- alpha_vec[m]
    
    if (beta_val > 1e-10) {
      if (beta_val < alp) {
        e1 <- (gamma(1 + 1 / alp) * gamma(1 / beta_val - 1 / alp)) / (lam * (beta_val^(1 / alp)) * gamma(1 / beta_val))
      } else {
        e1 <- Inf
      }
      
      if (beta_val < alp / 2) {
        e2 <- (gamma(1 + 2 / alp) * gamma(1 / beta_val - 2 / alp)) / ((lam^2) * (beta_val^(2 / alp)) * gamma(1 / beta_val))
        v <- e2 - e1^2
      } else {
        v <- Inf
      }
    } else {
      e1 <- gamma(1 + 1 / alp) / lam
      e2 <- gamma(1 + 2 / alp) / (lam^2)
      v <- e2 - e1^2
    }
    
    mean_vec[m] <- e1
    var_vec[m] <- v
  }
  
  data.frame(
    component = object$components,
    lambda = lambda_vec,
    alpha = alpha_vec,
    beta = rep(beta_val, M),
    mean = mean_vec,
    variance = var_vec,
    sd = sqrt(pmax(var_vec, 0))
  )
}
