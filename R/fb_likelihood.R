#' Log-Likelihood Evaluation for Frailty Models
#'
#' Computes the log-likelihood and individual log-likelihood contributions for
#' shared multiplicative, shared additive, univariate multiplicative, univariate additive,
#' and merged frailty models.
#'
#' @param params Numeric vector of all model parameters (baseline, regression coefficients, frailty).
#' @param data A list containing processed data elements (times, status, covariates).
#' @param baseline Baseline hazard distribution name.
#' @param frailty Frailty distribution name.
#' @param type Model type: \code{"multiplicative_bivariate"}, \code{"additive_bivariate"},
#'   \code{"multiplicative_univariate"}, \code{"additive_univariate"}, or \code{"merged_univariate"}.
#' @param return_indiv Logical; if \code{TRUE}, returns a vector of individual log-likelihoods.
#'
#' @return Total log-likelihood scalar or vector of individual log-likelihoods.
#' @export
#' @keywords internal
loglik_frailty <- function(params, data, baseline = "weibull", frailty = "gamma",
                           type = "multiplicative_bivariate", return_indiv = FALSE) {
  eps <- 1e-15
  pen_val <- -1e10
  
  # Check parameter validity
  n_base <- get_n_base_params(baseline)
  n_frail <- get_n_frail_params(frailty)
  n_cov <- data$p_cov
  
  if (length(params) != (n_base + n_cov + n_frail)) {
    if (return_indiv) return(rep(pen_val, data$n))
    return(pen_val)
  }
  
  base_params <- params[1:n_base]
  beta <- if (n_cov > 0) params[(n_base + 1):(n_base + n_cov)] else numeric(0)
  frail_params <- params[(n_base + n_cov + 1):(n_base + n_cov + n_frail)]
  
  # Ensure strictly positive parameters where required
  if (any(base_params <= 0) || any(frail_params <= 0)) {
    if (return_indiv) return(rep(pen_val, data$n))
    return(pen_val)
  }
  if (frailty == "posstable" && (frail_params[1] <= 0 || frail_params[1] >= 1)) {
    if (return_indiv) return(rep(pen_val, data$n))
    return(pen_val)
  }
  if (frailty == "mixture_invgauss" && (frail_params[1] <= 0 || frail_params[1] >= 1)) {
    if (return_indiv) return(rep(pen_val, data$n))
    return(pen_val)
  }
  
  log_lik_i <- numeric(data$n)
  
  if (type == "multiplicative_bivariate") {
    t1 <- data$time1
    t2 <- data$time2
    d1 <- data$status1
    d2 <- data$status2
    
    b1 <- eval_baseline(t1, baseline, base_params)
    b2 <- eval_baseline(t2, baseline, base_params)
    
    lp1 <- if (n_cov > 0) as.vector(data$X1 %*% beta) else rep(0, data$n)
    lp2 <- if (n_cov > 0) as.vector(data$X2 %*% beta) else rep(0, data$n)
    lp1 <- pmin(pmax(lp1, -30), 30)
    lp2 <- pmin(pmax(lp2, -30), 30)
    rho1 <- exp(lp1)
    rho2 <- exp(lp2)
    
    eta1 <- b1$H0 * rho1
    eta2 <- b2$H0 * rho2
    s <- eta1 + eta2
    
    fr <- eval_laplace(s, frailty, frail_params)
    
    c00 <- (d1 == 0 & d2 == 0)
    c10 <- (d1 == 1 & d2 == 0)
    c01 <- (d1 == 0 & d2 == 1)
    c11 <- (d1 == 1 & d2 == 1)
    
    if (any(c00)) log_lik_i[c00] <- fr$log_L0[c00]
    if (any(c10)) log_lik_i[c10] <- b1$log_h0[c10] + lp1[c10] + fr$log_L1[c10]
    if (any(c01)) log_lik_i[c01] <- b2$log_h0[c01] + lp2[c01] + fr$log_L1[c01]
    if (any(c11)) log_lik_i[c11] <- b1$log_h0[c11] + lp1[c11] + b2$log_h0[c11] + lp2[c11] + fr$log_L2[c11]
    
  } else if (type == "additive_bivariate") {
    t1 <- data$time1
    t2 <- data$time2
    d1 <- data$status1
    d2 <- data$status2
    
    b1 <- eval_baseline(t1, baseline, base_params)
    b2 <- eval_baseline(t2, baseline, base_params)
    
    lp1 <- if (n_cov > 0) as.vector(data$X1 %*% beta) else rep(0, data$n)
    lp2 <- if (n_cov > 0) as.vector(data$X2 %*% beta) else rep(0, data$n)
    lp1 <- pmin(pmax(lp1, -30), 30)
    lp2 <- pmin(pmax(lp2, -30), 30)
    rho1 <- exp(lp1)
    rho2 <- exp(lp2)
    
    zeta1 <- t1 * rho1
    zeta2 <- t2 * rho2
    s <- zeta1 + zeta2
    
    fr <- eval_laplace(s, frailty, frail_params)
    base_surv_log <- -(b1$H0 + b2$H0)
    
    c00 <- (d1 == 0 & d2 == 0)
    c10 <- (d1 == 1 & d2 == 0)
    c01 <- (d1 == 0 & d2 == 1)
    c11 <- (d1 == 1 & d2 == 1)
    
    if (any(c00)) {
      log_lik_i[c00] <- base_surv_log[c00] + fr$log_L0[c00]
    }
    if (any(c10)) {
      term10 <- pmax(b1$h0[c10] * fr$L0[c10] + rho1[c10] * fr$L1[c10], eps)
      log_lik_i[c10] <- base_surv_log[c10] + log(term10)
    }
    if (any(c01)) {
      term01 <- pmax(b2$h0[c01] * fr$L0[c01] + rho2[c01] * fr$L1[c01], eps)
      log_lik_i[c01] <- base_surv_log[c01] + log(term01)
    }
    if (any(c11)) {
      term11 <- pmax(b1$h0[c11] * b2$h0[c11] * fr$L0[c11] + 
                     (b1$h0[c11] * rho2[c11] + b2$h0[c11] * rho1[c11]) * fr$L1[c11] + 
                     rho1[c11] * rho2[c11] * fr$L2[c11], eps)
      log_lik_i[c11] <- base_surv_log[c11] + log(term11)
    }
    
  } else if (type == "multiplicative_univariate") {
    t <- data$time
    d <- data$status
    b <- eval_baseline(t, baseline, base_params)
    lp <- if (n_cov > 0) as.vector(data$X %*% beta) else rep(0, data$n)
    lp <- pmin(pmax(lp, -30), 30)
    rho <- exp(lp)
    
    eta <- b$H0 * rho
    fr <- eval_laplace(eta, frailty, frail_params)
    
    c0 <- (d == 0)
    c1 <- (d == 1)
    if (any(c0)) log_lik_i[c0] <- fr$log_L0[c0]
    if (any(c1)) log_lik_i[c1] <- b$log_h0[c1] + lp[c1] + fr$log_L1[c1]
    
  } else if (type == "additive_univariate") {
    t <- data$time
    d <- data$status
    b <- eval_baseline(t, baseline, base_params)
    lp <- if (n_cov > 0) as.vector(data$X %*% beta) else rep(0, data$n)
    lp <- pmin(pmax(lp, -30), 30)
    rho <- exp(lp)
    
    zeta <- t * rho
    fr <- eval_laplace(zeta, frailty, frail_params)
    
    c0 <- (d == 0)
    c1 <- (d == 1)
    if (any(c0)) log_lik_i[c0] <- -b$H0[c0] + fr$log_L0[c0]
    if (any(c1)) {
      term1 <- pmax(b$h0[c1] * fr$L0[c1] + rho[c1] * fr$L1[c1], eps)
      log_lik_i[c1] <- -b$H0[c1] + log(term1)
    }
  } else if (type == "merged_univariate") {
    # Merged multiplicative and additive frailty model
    # h(t | z1, z2) = [h0(t) + z1 * rho] * z2
    # S(t | z1) = L_Z2(H0(t) + z1 * t * rho)
    # S(t) = E_Z1[ L_Z2(H0(t) + Z1 * t * rho) ]
    # Solved using Gauss-Laguerre quadrature (20 nodes)
    t <- data$time
    d <- data$status
    b <- eval_baseline(t, baseline, base_params)
    lp <- if (n_cov > 0) as.vector(data$X %*% beta) else rep(0, data$n)
    lp <- pmin(pmax(lp, -30), 30)
    rho <- exp(lp)
    
    # Gauss-Laguerre 15 nodes and weights
    gl_nodes <- c(0.0933078, 0.4926917, 1.2155954, 2.2699495, 3.6676227,
                  5.4253366, 7.5659162, 10.1202286, 13.1302825, 16.6544077,
                  20.7764789, 25.6238942, 31.4075192, 38.5306833, 48.0260856)
    gl_weights <- c(0.2373016, 0.2862016, 0.2262832, 0.1333060, 0.0583818,
                    0.0187809, 0.0042724, 0.0006728, 0.0000703, 0.0000045,
                    0.00000017, 3.4e-9, 3.0e-11, 8.8e-14, 3.1e-17)
    
    theta <- frail_params[1]
    # S(t_i) and f(t_i)
    S_i <- numeric(data$n)
    f_i <- numeric(data$n)
    
    for (k in seq_along(gl_nodes)) {
      # z1 = gl_nodes[k] / theta (since z1 ~ Gamma(theta, theta))
      z1_val <- gl_nodes[k] / theta
      # node weight adjustment for Gamma(theta, theta): (theta^theta / Gamma(theta)) * z1^(theta-1) * exp(-theta z1)
      # with substitution u = theta * z1 = gl_nodes[k], du = theta dz1,
      # integral over z1 becomes integral over u of u^(theta-1) / Gamma(theta) * exp(-u)
      wk_adj <- gl_weights[k] * (gl_nodes[k]^(theta - 1)) / gamma(theta)
      if (is.finite(wk_adj) && wk_adj > 0) {
        s_val <- b$H0 + z1_val * t * rho
        # L_Z2 and derivatives
        u2 <- 1 + s_val / theta
        surv_cond <- u2^(-theta)
        haz_multiplier <- (b$h0 + z1_val * rho) * (u2^(-theta - 1))
        
        S_i <- S_i + wk_adj * surv_cond
        f_i <- f_i + wk_adj * haz_multiplier
      }
    }
    
    S_i <- pmin(pmax(S_i, eps), 1)
    f_i <- pmax(f_i, eps)
    
    log_lik_i <- ifelse(d == 1, log(f_i), log(S_i))
  }
  
  # Numerical safety
  log_lik_i[!is.finite(log_lik_i)] <- pen_val
  log_lik_i <- pmax(log_lik_i, pen_val)
  
  if (return_indiv) return(log_lik_i)
  sum(log_lik_i)
}

#' Number of Baseline Parameters
#' @param baseline Baseline name.
#' @return Integer count.
#' @keywords internal
get_n_base_params <- function(baseline) {
  if (baseline == "exponential") return(1)
  if (baseline %in% c("weibull", "genexponential", "exppower", "lfr", 
                      "genrayleigh", "weightedexp", "extendedweibull")) return(2)
  if (baseline %in% c("genloglogistic", "genweibull")) return(3)
  2
}

#' Number of Frailty Parameters
#' @param frailty Frailty name.
#' @return Integer count.
#' @keywords internal
get_n_frail_params <- function(frailty) {
  if (frailty %in% c("gamma", "invgauss", "posstable")) return(1)
  if (frailty %in% c("genlindley", "weightedlindley")) return(2)
  if (frailty == "compoundpoisson") return(2)
  if (frailty == "compoundnb") return(3)
  if (frailty == "mixture_invgauss") return(3)
  1
}
