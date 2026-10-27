#' Simulate Dependent Component Lifetimes via Clayton Copula
#'
#' @description
#' Implements the conditional sampling algorithm using the Clayton survival copula
#' (Equations 16--17 of Ling et al., 2021) to generate correlated component failure
#' times under Lomax distributions.
#'
#' @param n An integer specifying the number of devices to generate.
#' @param lambda A positive numeric vector of failure rates \eqn{\lambda_m} for each component (\eqn{M = \text{length}(\lambda)}).
#' @param beta A positive numeric scalar representing the frailty variance / copula dependence parameter (\eqn{\beta > 0}).
#' @param seed An optional integer for reproducibility.
#'
#' @return An \eqn{n \times M} numeric matrix of simulated component failure times,
#'   where each row represents a device and each column represents a component.
#'
#' @details
#' The algorithm generates \eqn{u_1, v_2, \dots, v_M \sim U(0, 1)} and successively computes:
#' \deqn{u_2 = \left( 1 - (1 - v_2)^{-\frac{\beta}{\beta+1}} u_1^{-\beta} \right)^{-1/\beta}}
#' \deqn{u_m = \left( \left( v_m^{-\frac{\beta}{(m-1)\beta+1}} - 1 \right) \sum_{j=1}^{m-1} u_j^{-\beta} - (m - 2) v_m^{-\frac{\beta}{(m-1)\beta+1}} + (m - 1) \right)^{-1/\beta}}
#' and component failure times \eqn{t_m = \frac{u_m^{-\beta} - 1}{\beta \lambda_m}}.
#'
#' Missing values are not permitted.
#'
#' @references
#' Ling, M.-H., Balakrishnan, N., Yu, C., and So, H. Y. (2021). Inference for One-Shot
#' Devices with Dependent k-Out-of-M Structured Components under Gamma Frailty.
#' \emph{Mathematics}, 9(23), 3032. \doi{10.3390/math9233032}.
#'
#' @examples
#' lam <- c(0.0086, 0.0067, 0.0052, 0.0024)
#' sim_times <- sim_copula(n = 100, lambda = lam, beta = 0.3, seed = 123)
#' head(sim_times)
#' @export
sim_copula <- function(n, lambda, beta, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  if (missing(n) || is.null(n) || length(n) != 1 || is.na(n) || n < 1) {
    stop("'n' must be a positive integer.")
  }
  if (missing(lambda) || is.null(lambda) || any(is.na(lambda)) || any(lambda <= 0)) {
    stop("'lambda' must be a positive numeric vector.")
  }
  if (missing(beta) || is.null(beta) || length(beta) != 1 || is.na(beta) || beta <= 0) {
    stop("'beta' must be a positive numeric scalar.")
  }
  
  M <- length(lambda)
  n <- as.integer(n)
  
  U_mat <- matrix(0, nrow = n, ncol = M)
  U_mat[, 1] <- stats::runif(n)
  
  if (M >= 2) {
    v2 <- stats::runif(n)
    term2 <- 1 - (1 - v2)^(-beta / (beta + 1)) * (U_mat[, 1]^(-beta))
    term2[term2 <= 0] <- 1e-10
    U_mat[, 2] <- term2^(-1 / beta)
  }
  
  if (M >= 3) {
    for (m in 3:M) {
      vm <- stats::runif(n)
      gamma_pow <- -beta / ((m - 1) * beta + 1)
      vm_pow <- vm^gamma_pow
      
      sum_u_beta <- rowSums(U_mat[, 1:(m - 1), drop = FALSE]^(-beta))
      term_m <- (vm_pow - 1) * sum_u_beta - (m - 2) * vm_pow + (m - 1)
      term_m[term_m <= 0] <- 1e-10
      U_mat[, m] <- term_m^(-1 / beta)
    }
  }
  
  T_mat <- matrix(0, nrow = n, ncol = M)
  for (m in seq_len(M)) {
    T_mat[, m] <- (U_mat[, m]^(-beta) - 1) / (beta * lambda[m])
    T_mat[, m] <- pmax(T_mat[, m], 0)
  }
  
  colnames(T_mat) <- paste0("T_", seq_len(M))
  T_mat
}

#' Simulate Dependent Component Lifetimes via Gamma Frailty
#'
#' @description
#' Generates correlated component lifetimes by first generating a latent gamma frailty
#' \eqn{\gamma \sim \text{Gamma}(1/\beta, \beta)} and then simulating conditionally
#' independent Exponential or Weibull component lifetimes.
#'
#' @param n An integer specifying the number of devices to generate.
#' @param lambda A positive numeric vector of failure rates \eqn{\lambda_m} for each component (\eqn{M = \text{length}(\lambda)}).
#' @param beta A positive numeric scalar representing the frailty variance (\eqn{\beta > 0}).
#' @param alpha An optional positive numeric vector of Weibull shape parameters (default 1 for exponential).
#' @param seed An optional integer for reproducibility.
#'
#' @return An \eqn{n \times M} numeric matrix of simulated component failure times.
#'
#' @details
#' Conditioned on \eqn{\gamma_{i,j} \sim \text{Gamma}(\text{shape} = 1/\beta, \text{scale} = \beta)},
#' the component lifetimes are independent with:
#' \deqn{T_{i,j,m} = \frac{1}{\lambda_m} \left( -\frac{\ln(1 - u_m)}{\gamma_{i,j}} \right)^{1/\alpha_m}}
#' where \eqn{u_m \sim U(0, 1)}.
#'
#' Missing values are not permitted.
#'
#' @references
#' Ling, M.-H., Balakrishnan, N., Yu, C., and So, H. Y. (2021). Inference for One-Shot
#' Devices with Dependent k-Out-of-M Structured Components under Gamma Frailty.
#' \emph{Mathematics}, 9(23), 3032. \doi{10.3390/math9233032}.
#'
#' @examples
#' lam <- c(0.0086, 0.0067, 0.0052, 0.0024)
#' sim_times <- sim_frailty(n = 100, lambda = lam, beta = 0.3, seed = 123)
#' head(sim_times)
#' @export
sim_oneshot_frailty_lifetimes <- function(n, lambda, beta, alpha = NULL, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  if (missing(n) || is.null(n) || length(n) != 1 || is.na(n) || n < 1) {
    stop("'n' must be a positive integer.")
  }
  if (missing(lambda) || is.null(lambda) || any(is.na(lambda)) || any(lambda <= 0)) {
    stop("'lambda' must be a positive numeric vector.")
  }
  if (missing(beta) || is.null(beta) || length(beta) != 1 || is.na(beta) || beta <= 0) {
    stop("'beta' must be a positive numeric scalar.")
  }
  
  M <- length(lambda)
  n <- as.integer(n)
  if (is.null(alpha)) alpha <- rep(1, M)
  
  # Latent frailty: shape = 1/beta, scale = beta (mean = 1, var = beta)
  gamma_vec <- stats::rgamma(n, shape = 1 / beta, scale = beta)
  
  T_mat <- matrix(0, nrow = n, ncol = M)
  for (m in seq_len(M)) {
    u_m <- stats::runif(n)
    T_mat[, m] <- (1 / lambda[m]) * (-log(1 - u_m) / gamma_vec)^(1 / alpha[m])
  }
  
  colnames(T_mat) <- paste0("T_", seq_len(M))
  T_mat
}

#' Simulate One-Shot Device Testing Datasets
#'
#' @description
#' Generates synthetic one-shot device testing data frames across multiple test groups
#' under constant-stress accelerated life-testing plans with dependent components.
#'
#' @param n_groups An integer specifying the number of test groups \eqn{I}.
#' @param K An integer scalar or vector of length \eqn{I} specifying sample sizes per group.
#' @param inspection_times A positive numeric vector of inspection times \eqn{\tau_i}.
#' @param stress_levels A numeric vector or matrix of stress levels \eqn{s_i}.
#' @param a_params An \eqn{M \times p} matrix of scale regression coefficients.
#' @param b_params An optional \eqn{M \times p} matrix of shape regression coefficients (for Weibull).
#' @param beta A positive numeric scalar representing the frailty variance (\eqn{\beta > 0}).
#' @param baseline A character string: \code{"exponential"} (default) or \code{"weibull"}.
#' @param method A character string: \code{"frailty"} (default) or \code{"copula"}.
#' @param seed An optional integer for reproducibility.
#'
#' @return A data frame containing:
#'   \describe{
#'     \item{group}{Test group index \eqn{i}.}
#'     \item{stress}{Stress level \eqn{s_i}.}
#'     \item{time}{Inspection time \eqn{\tau_i}.}
#'     \item{counts}{Observed counts for each of the \eqn{2^M} malfunction combinations.}
#'   }
#'
#' @details
#' Missing values are not permitted in any arguments.
#'
#' @references
#' Ling, M.-H., Balakrishnan, N., Yu, C., and So, H. Y. (2021). Inference for One-Shot
#' Devices with Dependent k-Out-of-M Structured Components under Gamma Frailty.
#' \emph{Mathematics}, 9(23), 3032. \doi{10.3390/math9233032}.
#'
#' @examples
#' a_par <- matrix(c(-6.0, 0.05, -6.5, 0.06), nrow = 2, byrow = TRUE)
#' sim_df <- sim_oneshot(
#'   n_groups = 4,
#'   K = 50,
#'   inspection_times = c(10, 20, 10, 20),
#'   stress_levels = c(35, 35, 45, 45),
#'   a_params = a_par,
#'   beta = 0.2,
#'   seed = 42
#' )
#' head(sim_df)
#' @export
sim_oneshot <- function(n_groups = 4, K = 50, inspection_times = c(1000, 2000, 3000, 4000), stress_levels = c(25, 35, 45, 55), a_params,
                        b_params = NULL, beta = 0.2,
                        baseline = c("exponential", "weibull"),
                        method = c("frailty", "copula"),
                        seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  baseline <- match.arg(baseline)
  method <- match.arg(method)
  
  I_groups <- as.integer(n_groups)
  if (length(K) == 1) K <- rep(as.integer(K), I_groups)
  if (length(inspection_times) == 1) inspection_times <- rep(inspection_times, I_groups)
  
  if (is.vector(stress_levels)) {
    stress_mat <- matrix(stress_levels, ncol = 1)
  } else {
    stress_mat <- as.matrix(stress_levels)
  }
  
  X_cov <- cbind(1, stress_mat)
  p_cov <- ncol(X_cov)
  M <- nrow(a_params)
  powerset <- generate_powerset(M)
  n_subsets <- length(powerset)
  
  col_names <- vapply(powerset, function(s) {
    if (length(s) == 0) "n_none"
    else if (length(s) == M) "n_all"
    else paste0("n_", paste(s, collapse = "_"))
  }, character(1))
  
  count_matrix <- matrix(0, nrow = I_groups, ncol = n_subsets)
  
  for (i in seq_len(I_groups)) {
    x_i <- matrix(X_cov[i, ], nrow = 1)
    lam_vec <- numeric(M)
    alpha_vec <- numeric(M)
    for (m in seq_len(M)) {
      lam_vec[m] <- exp(sum(x_i * a_params[m, ]))
      alpha_vec[m] <- if (baseline == "weibull" && !is.null(b_params)) exp(sum(x_i * b_params[m, ])) else 1
    }
    
    # Simulate lifetimes for K[i] devices
    if (method == "frailty") {
      T_sample <- sim_oneshot_frailty_lifetimes(n = K[i], lambda = lam_vec, beta = beta, alpha = alpha_vec)
    } else {
      T_sample <- sim_copula(n = K[i], lambda = lam_vec, beta = beta)
    }
    
    # Determine malfunctioned subset for each device at inspection time tau_i
    tau_i <- inspection_times[i]
    failed_mat <- (T_sample <= tau_i) # TRUE if failed before tau_i
    
    for (d_idx in seq_len(K[i])) {
      failed_comps <- which(failed_mat[d_idx, ])
      
      # Match with powerset
      matched_idx <- 1
      for (s_idx in seq_along(powerset)) {
        if (length(powerset[[s_idx]]) == length(failed_comps) &&
            all(sort(powerset[[s_idx]]) == sort(failed_comps))) {
          matched_idx <- s_idx
          break
        }
      }
      count_matrix[i, matched_idx] <- count_matrix[i, matched_idx] + 1
    }
  }
  
  df_out <- data.frame(
    group = seq_len(I_groups),
    stress = stress_mat[, 1],
    time = inspection_times,
    stringsAsFactors = FALSE
  )
  
  count_df <- as.data.frame(count_matrix)
  colnames(count_df) <- col_names
  
  res_df <- cbind(df_out, count_df)
  if ("n_all" %in% colnames(res_df) && !"n_both" %in% colnames(res_df)) {
    res_df$n_both <- res_df$n_all
  }
  res_df
}
