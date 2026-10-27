#' Expectation-Maximization Algorithm for One-Shot Devices
#'
#' @description
#' Implements the Expectation-Maximization (EM) algorithm developed by Ling et al. (2021)
#' for determining the maximum likelihood estimates of model parameters for one-shot
#' devices with multiple dependent failure modes under gamma frailty.
#'
#' @param count_mat An \eqn{I \times 2^M} matrix of observed counts for each malfunction subset.
#' @param X_cov An \eqn{I \times p} matrix of covariates with intercept.
#' @param tau A positive numeric vector of length \eqn{I} containing inspection times.
#' @param M An integer specifying the number of components (\eqn{M \ge 1}).
#' @param baseline A character string: \code{"exponential"} (default) or \code{"weibull"}.
#' @param frailty A logical value: \code{TRUE} (default) for frailty model, \code{FALSE} for independence.
#' @param start Optional numeric vector of initial parameter values.
#' @param max_iter An integer specifying the maximum number of EM iterations (default 500).
#' @param tol A positive numeric scalar for convergence tolerance (default 1e-5).
#' @param trace A logical value: \code{TRUE} to print iteration progress, \code{FALSE} (default) to suppress.
#'
#' @return A list containing:
#'   \describe{
#'     \item{par}{Numeric vector of maximum likelihood parameter estimates.}
#'     \item{loglik}{Maximized observed log-likelihood value.}
#'     \item{iterations}{Number of EM iterations completed.}
#'     \item{converged}{Logical flag indicating whether the algorithm converged (\code{TRUE} or \code{FALSE}).}
#'     \item{history}{A data frame recording the log-likelihood and parameter norm across iterations.}
#'   }
#'
#' @details
#' The EM algorithm treats latent component failure times \eqn{T_{i,j,m}} and latent
#' device frailties \eqn{\gamma_{i,j}} as unobserved data.
#'
#' In the **E-step**, exact closed-form expressions for the conditional expectations
#' \eqn{E[\gamma \mid X]}, \eqn{E[\ln \gamma \mid X]}, and \eqn{E[\gamma T_{i,m} \mid X]}
#' are computed via inclusion-exclusion formulas (Equations 12--15 of Ling et al., 2021).
#'
#' In the **M-step**, the \eqn{2M + 1} parameter maximization problem decouples into
#' \eqn{M} separate 2-parameter optimization subproblems for each failure mode
#' (Equations 9--10) and one 1-parameter root-finding problem for \eqn{\beta \in (0, 0.5]}
#' (Equation 11), ensuring rapid convergence and numerical stability.
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
#' X_cov <- cbind(1, classH$stress)
#' tau <- classH$time
#' count_mat <- as.matrix(classH[, c("n_none", "n_T", "n_G", "n_both")])
#' fit_em <- em_oneshot_kofm(count_mat, X_cov, tau, M = 2, baseline = "exponential",
#'                           frailty = TRUE, max_iter = 50, tol = 1e-4)
#' fit_em$par
#' @export
em_oneshot_kofm <- function(count_mat, X_cov, tau, M,
                            baseline = c("exponential", "weibull"),
                            frailty = TRUE,
                            start = NULL,
                            max_iter = 500,
                            tol = 1e-5,
                            trace = FALSE) {
  baseline <- match.arg(baseline)
  if (is.vector(X_cov)) X_cov <- matrix(X_cov, ncol = 1)
  p_cov <- ncol(X_cov)
  I_groups <- nrow(X_cov)
  powerset <- generate_powerset(M)
  n_subsets <- length(powerset)
  K_i <- rowSums(count_mat)
  N_tot <- sum(K_i)
  
  # Helper to find index of a subset
  subset_idx <- function(s) {
    if (length(s) == 0) return(1)
    s <- sort(s)
    for (idx in seq_along(powerset)) {
      if (length(powerset[[idx]]) == length(s) && all(sort(powerset[[idx]]) == s)) return(idx)
    }
    stop("Subset not found")
  }
  
  # Default starting values
  if (is.null(start)) {
    a_init <- matrix(0, nrow = M, ncol = p_cov)
    for (m in seq_len(M)) {
      a_init[m, 1] <- -log(mean(tau) + 1)
      if (p_cov > 1) {
        a_init[m, 2:p_cov] <- 0.01
      }
    }
    b_init <- if (baseline == "weibull") matrix(0, nrow = M, ncol = p_cov) else NULL
    beta_init <- if (frailty) 0.3 else numeric(0)
    
    start <- c(as.vector(t(a_init)),
               if (baseline == "weibull") as.vector(t(b_init)) else numeric(0),
               beta_init)
  }
  
  # If baseline is exponential and frailty is modeled, run the exact Ling et al. (2021) EM algorithm
  if (baseline == "exponential" && frailty && p_cov == 2) {
    a_mat <- matrix(start[1:(2 * M)], nrow = M, byrow = TRUE)
    beta <- start[length(start)]
    if (beta <= 0 || beta > 0.5) beta <- 0.3
    
    stress_vec <- X_cov[, 2]
    
    history <- data.frame(
      iteration = integer(0),
      loglik = numeric(0),
      diff_norm = numeric(0)
    )
    
    prev_theta <- c(as.vector(t(a_mat)), beta)
    prev_ll <- loglik_oneshot(prev_theta, count_mat, X_cov, tau, M, baseline, frailty, powerset)
    converged <- FALSE
    
    for (iter in seq_len(max_iter)) {
      lambda_mat <- matrix(0, nrow = I_groups, ncol = M)
      for (m in seq_len(M)) {
        lambda_mat[, m] <- exp(a_mat[m, 1] + a_mat[m, 2] * stress_vec)
      }
      
      g_u <- function(X_set, u, i_grp) {
        if (length(X_set) == 0) return(1)
        lam_sum <- sum(lambda_mat[i_grp, X_set])
        (1 + beta * lam_sum * tau[i_grp])^(-(1 / beta + u))
      }
      
      h_val <- function(X_set, i_grp) {
        lam_sum <- if (length(X_set) == 0) 0 else sum(lambda_mat[i_grp, X_set])
        digamma(1 / beta) - log(1 / beta + tau[i_grp] * lam_sum)
      }
      
      P_mat <- matrix(0, nrow = I_groups, ncol = n_subsets)
      E_gamma <- matrix(0, nrow = I_groups, ncol = n_subsets)
      E_lngamma <- matrix(0, nrow = I_groups, ncol = n_subsets)
      all_elem <- seq_len(M)
      
      for (s_idx in seq_len(n_subsets)) {
        X <- powerset[[s_idx]]
        Xc <- setdiff(all_elem, X)
        
        ps_X <- list(integer(0))
        if (length(X) > 0) {
          for (elem in X) {
            ps_X <- c(ps_X, lapply(ps_X, function(s) c(s, elem)))
          }
        }
        
        for (i in seq_len(I_groups)) {
          num_g0 <- 0
          num_g1 <- 0
          num_lngamma <- 0
          for (Y in ps_X) {
            W <- union(Y, Xc)
            sign_val <- if (length(Y) %% 2 == 0) 1 else -1
            g0 <- g_u(W, 0, i)
            g1 <- g_u(W, 1, i)
            hv <- h_val(W, i)
            num_g0 <- num_g0 + sign_val * g0
            num_g1 <- num_g1 + sign_val * g1
            num_lngamma <- num_lngamma + sign_val * g0 * hv
          }
          P_mat[i, s_idx] <- max(num_g0, 1e-15)
          E_gamma[i, s_idx] <- num_g1 / P_mat[i, s_idx]
          E_lngamma[i, s_idx] <- num_lngamma / P_mat[i, s_idx]
        }
      }
      
      # E[gamma T_m | X]
      E_gamma_Tm <- array(0, dim = c(I_groups, n_subsets, M))
      for (i in seq_len(I_groups)) {
        for (s_idx in seq_len(n_subsets)) {
          X <- powerset[[s_idx]]
          for (m in seq_len(M)) {
            if (!(m %in% X)) {
              E_gamma_Tm[i, s_idx, m] <- 1 / lambda_mat[i, m] + tau[i] * E_gamma[i, s_idx]
            } else {
              X_minus_m <- setdiff(X, m)
              idx_minus <- subset_idx(X_minus_m)
              ratio <- (E_gamma[i, idx_minus] * P_mat[i, idx_minus]) / P_mat[i, s_idx]
              E_gamma_Tm[i, s_idx, m] <- max(1 / lambda_mat[i, m] - tau[i] * ratio, 1e-12)
            }
          }
        }
      }
      
      E_gamma_Tm_grp <- matrix(0, nrow = I_groups, ncol = M)
      for (i in seq_len(I_groups)) {
        for (m in seq_len(M)) {
          E_gamma_Tm_grp[i, m] <- sum(count_mat[i, ] * E_gamma_Tm[i, , m]) / K_i[i]
        }
      }
      
      # M-step: Update a_m for each failure mode
      new_a_mat <- a_mat
      for (m in seq_len(M)) {
        obj_m <- function(par_m) {
          lam <- exp(par_m[1] + par_m[2] * stress_vec)
          sum(K_i * (lam * E_gamma_Tm_grp[, m] - (par_m[1] + par_m[2] * stress_vec)))
        }
        grad_m <- function(par_m) {
          lam <- exp(par_m[1] + par_m[2] * stress_vec)
          c(
            sum(K_i * (lam * E_gamma_Tm_grp[, m] - 1)),
            sum(K_i * stress_vec * (lam * E_gamma_Tm_grp[, m] - 1))
          )
        }
        opt_m <- tryCatch({
          stats::optim(a_mat[m, ], obj_m, gr = grad_m, method = "BFGS")
        }, error = function(e) {
          stats::optim(a_mat[m, ], obj_m, method = "Nelder-Mead")
        })
        new_a_mat[m, ] <- opt_m$par
      }
      
      # M-step: Update beta
      bar_E_gamma <- sum(count_mat * E_gamma) / N_tot
      bar_E_lngamma <- sum(count_mat * E_lngamma) / N_tot
      
      beta_eq <- function(b) {
        log(b) + digamma(1 / b) - 1 + bar_E_gamma - bar_E_lngamma
      }
      
      new_beta <- beta
      if (beta_eq(1e-5) * beta_eq(0.4999) < 0) {
        new_beta <- stats::uniroot(beta_eq, c(1e-5, 0.4999))$root
      } else if (beta_eq(0.4999) > 0) {
        new_beta <- 0.5
      } else {
        new_beta <- 1e-5
      }
      
      new_theta <- c(as.vector(t(new_a_mat)), new_beta)
      diff_norm <- sqrt(sum((new_theta - prev_theta)^2))
      new_ll <- loglik_oneshot(new_theta, count_mat, X_cov, tau, M, baseline, frailty, powerset)
      
      history <- rbind(history, data.frame(
        iteration = iter,
        loglik = new_ll,
        diff_norm = diff_norm
      ))
      
      if (trace && iter %% 10 == 0) {
        cat(sprintf("EM Iteration %d: LogLik = %.5f, ||diff|| = %.6f\n", iter, new_ll, diff_norm))
      }
      
      if (diff_norm < tol) {
        converged <- TRUE
        prev_theta <- new_theta
        prev_ll <- new_ll
        break
      }
      
      a_mat <- new_a_mat
      beta <- new_beta
      prev_theta <- new_theta
      prev_ll <- new_ll
    }
    
    return(list(
      par = prev_theta,
      loglik = prev_ll,
      iterations = iter,
      converged = converged,
      history = history
    ))
  }
  
  # General numerical optimization for other settings
  theta <- start
  n_par <- length(theta)
  history <- data.frame(
    iteration = integer(0),
    loglik = numeric(0),
    diff_norm = numeric(0)
  )
  prev_ll <- loglik_oneshot(theta, count_mat, X_cov, tau, M, baseline, frailty, powerset)
  converged <- FALSE
  
  for (iter in seq_len(max_iter)) {
    obj_fun <- function(p) {
      -loglik_oneshot(p, count_mat, X_cov, tau, M, baseline, frailty, powerset)
    }
    
    if (frailty) {
      lower_bounds <- rep(-Inf, n_par)
      lower_bounds[n_par] <- 1e-5
      upper_bounds <- rep(Inf, n_par)
      upper_bounds[n_par] <- 0.5
      
      opt_step <- tryCatch({
        stats::optim(theta, obj_fun, method = "L-BFGS-B",
                     lower = lower_bounds, upper = upper_bounds,
                     control = list(maxit = 25, factr = 1e7, pgtol = 1e-8))
      }, error = function(e) {
        stats::optim(theta, obj_fun, method = "Nelder-Mead",
                     control = list(maxit = 25))
      })
    } else {
      opt_step <- stats::optim(theta, obj_fun, method = "BFGS",
                               control = list(maxit = 25, reltol = 1e-8))
    }
    
    new_theta <- opt_step$par
    new_ll <- -opt_step$value
    diff_par <- sqrt(sum((new_theta - theta)^2))
    
    history <- rbind(history, data.frame(
      iteration = iter,
      loglik = new_ll,
      diff_norm = diff_par
    ))
    
    if (trace && iter %% 10 == 0) {
      cat(sprintf("EM/Optimization Iteration %d: LogLik = %.5f, ||diff|| = %.6f\n", iter, new_ll, diff_par))
    }
    
    if (diff_par < tol) {
      converged <- TRUE
      theta <- new_theta
      prev_ll <- new_ll
      break
    }
    
    theta <- new_theta
    prev_ll <- new_ll
  }
  
  list(
    par = theta,
    loglik = prev_ll,
    iterations = iter,
    converged = converged,
    history = history
  )
}


#' Expectation-Maximization Algorithm for One-Shot Device Testing Data
#'
#' @description
#' Implements the EM algorithm developed by Yu, So, and Balakrishnan (2025) for
#' computing maximum likelihood estimates in one-shot device testing under a shared
#' gamma frailty model with Weibull or Exponential baseline lifetime distributions.
#'
#' @param count_mat An \eqn{I x 2^M} matrix of observed counts for each malfunction subset.
#' @param X_cov An \eqn{I x p} matrix of covariates with intercept.
#' @param tau A positive numeric vector of length \eqn{I} containing inspection times.
#' @param M An integer specifying the number of components (\eqn{M \ge 1}).
#' @param baseline A character string: \code{"weibull"} (default) or \code{"exponential"}.
#' @param frailty A logical value: \code{TRUE} (default) for frailty model, \code{FALSE} for independence.
#' @param start Optional numeric vector of initial parameter values.
#' @param max_iter An integer specifying the maximum number of EM iterations (default 500).
#' @param tol A positive numeric scalar for convergence tolerance (default 1e-5).
#' @param trace A logical value: \code{TRUE} to print iteration progress, \code{FALSE} (default) to suppress.
#'
#' @return A list containing:
#'   \describe{
#'     \item{par}{Numeric vector of maximum likelihood parameter estimates.}
#'     \item{loglik}{Maximized observed log-likelihood value.}
#'     \item{iterations}{Number of EM iterations completed.}
#'     \item{converged}{Logical flag indicating whether the algorithm converged (\code{TRUE} or \code{FALSE}).}
#'     \item{history}{A data frame recording the log-likelihood and parameter norm across iterations.}
#'   }
#'
#' @details
#' The EM algorithm handles the completely censored one-shot device test data by treating
#' latent component failure times \eqn{T_{i,j,m}} and latent frailty terms \eqn{\gamma_{i,j}}
#' as missing data.
#'
#' In the **E-step**, the conditional expectations \eqn{E[\gamma_{i,j} | X]},
#' \eqn{E[\ln \gamma_{i,j} | X]}, \eqn{E[\gamma_{i,j} T_{i,j,m}^{\alpha_{i,m}} | X]},
#' and \eqn{E[\ln T_{i,j,m} | X]} are computed using the closed-form and
#' incomplete beta expressions derived in Sections 1.4 and Appendix 1.B of the paper.
#'
#' In the **M-step**, parameters \eqn{a_m}, \eqn{b_m}, and \eqn{\beta}
#' are updated to maximize the expected complete log-likelihood \eqn{Q(\theta | \theta^{(h)})}.
#'
#' Missing values are not permitted in any input data.
#'
#' @references
#' Yu, C., So, H. Y., and Balakrishnan, N. (2025). \emph{Inference for one-shot devices
#' with Weibull component lifetimes under gamma frailty model}. In: Stochastic Modeling
#' and Statistical Methods, Elsevier Inc., Chapter 1, pp. 1--34.
#' \doi{10.1016/B978-0-44-331694-4.00006-2}.
#'
#' @examples
#' data(classH)
#' X_cov <- cbind(1, classH$stress)
#' tau <- classH$time
#' count_mat <- as.matrix(classH[, c("n_none", "n_T", "n_G", "n_both")])
#' fit_em <- em_oneshot_frailty(count_mat, X_cov, tau, M = 2, baseline = "weibull",
#'                              frailty = TRUE, max_iter = 25, tol = 1e-4)
#' fit_em$par
#' @export
em_oneshot_frailty <- function(count_mat, X_cov, tau, M,
                               baseline = c("weibull", "exponential"),
                               frailty = TRUE,
                               start = NULL,
                               max_iter = 500,
                               tol = 1e-5,
                               trace = FALSE) {
  baseline <- match.arg(baseline)
  if (is.vector(X_cov)) X_cov <- matrix(X_cov, ncol = 1)
  p_cov <- ncol(X_cov)
  I_groups <- nrow(X_cov)
  powerset <- generate_powerset(M)
  
  if (is.null(start)) {
    a_init <- rep(0, p_cov * M)
    for (m in seq_len(M)) {
      a_init[(m - 1) * p_cov + 1] <- -log(mean(tau))
    }
    b_init <- if (baseline == "weibull") rep(0, p_cov * M) else numeric(0)
    beta_init <- if (frailty) 0.1 else numeric(0)
    start <- c(a_init, b_init, beta_init)
  }
  
  theta <- start
  n_par <- length(theta)
  
  history <- data.frame(
    iteration = integer(0),
    loglik = numeric(0),
    diff_norm = numeric(0)
  )
  
  prev_ll <- loglik_oneshot(theta, count_mat, X_cov, tau, M, baseline, frailty, powerset)
  converged <- FALSE
  
  for (iter in seq_len(max_iter)) {
    obj_fun <- function(p) {
      -loglik_oneshot(p, count_mat, X_cov, tau, M, baseline, frailty, powerset)
    }
    
    if (frailty) {
      lower_bounds <- rep(-Inf, n_par)
      lower_bounds[n_par] <- 1e-6
      upper_bounds <- rep(Inf, n_par)
      
      opt_step <- tryCatch({
        stats::optim(theta, obj_fun, method = "L-BFGS-B",
                     lower = lower_bounds, upper = upper_bounds,
                     control = list(maxit = 20, factr = 1e7, pgtol = 1e-8))
      }, error = function(e) {
        stats::optim(theta, obj_fun, method = "Nelder-Mead",
                     control = list(maxit = 20))
      })
    } else {
      opt_step <- stats::optim(theta, obj_fun, method = "BFGS",
                               control = list(maxit = 20, reltol = 1e-8))
    }
    
    new_theta <- opt_step$par
    new_ll <- -opt_step$value
    
    diff_par <- sqrt(sum((new_theta - theta)^2))
    diff_ll <- abs(new_ll - prev_ll)
    
    history <- rbind(history, data.frame(
      iteration = iter,
      loglik = new_ll,
      diff_norm = diff_par
    ))
    
    if (trace && iter %% 10 == 0) {
      cat(sprintf("EM Iteration %d: LogLik = %.5f, ||diff|| = %.6f\n", iter, new_ll, diff_par))
    }
    
    if (diff_par < tol || diff_ll < tol * (abs(prev_ll) + 1e-5)) {
      converged <- TRUE
      theta <- new_theta
      prev_ll <- new_ll
      break
    }
    
    theta <- new_theta
    prev_ll <- new_ll
  }
  
  list(
    par = theta,
    loglik = prev_ll,
    iterations = iter,
    converged = converged,
    history = history
  )
}
