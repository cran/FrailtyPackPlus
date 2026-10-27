#' Simulate Survival Data with Shared or Univariate Frailty
#'
#' Generates synthetic survival datasets under any supported baseline hazard
#' and frailty distribution with right censoring.
#'
#' @param n Sample size (number of subjects/clusters).
#' @param baseline Baseline hazard name (default \code{"weibull"}).
#' @param frailty Frailty distribution name (default \code{"gamma"}).
#' @param type Model type (default \code{"multiplicative_bivariate"}).
#' @param beta Regression coefficients vector.
#' @param theta Frailty parameter(s).
#' @param baseline_params Baseline parameters vector.
#' @param censor_rate Right-censoring probability.
#'
#' @return A data frame containing simulated survival times, status indicators,
#'   covariates, and true frailty values.
#' @export
#' @examples
#' sim_data <- sim_frailty(
#'   n = 50,
#'   baseline = "weibull",
#'   frailty = "gamma",
#'   type = "multiplicative_bivariate",
#'   beta = c(0.5, -0.5),
#'   theta = 1.5,
#'   baseline_params = c(0.02, 1.2)
#' )
#' head(sim_data)
fb_sim_frailty <- function(n = 100, baseline = "weibull", frailty = "gamma",
                        type = "multiplicative_bivariate", beta = c(0.5, -0.5),
                        theta = 1.0, baseline_params = c(0.01, 1.2),
                        censor_rate = 0.2) {
  z <- r_frailty(n, frailty, theta)
  x1 <- stats::rnorm(n, 0, 1)
  x2 <- stats::rbinom(n, 1, 0.5)
  X <- cbind(x1 = x1, x2 = x2)
  beta_vec <- rep(beta, length.out = ncol(X))
  lp <- as.vector(X %*% beta_vec)
  rho <- exp(lp)
  
  inv_H0 <- function(u, b_name, b_par) {
    if (b_name == "weibull") {
      (u / b_par[1])^(1 / b_par[2])
    } else if (b_name == "exponential") {
      u / b_par[1]
    } else if (b_name == "exppower") {
      (log(u + 1))^(1 / b_par[2]) / b_par[1]
    } else {
      sapply(u, function(u_val) {
        if (u_val <= 0) return(0)
        res <- tryCatch(
          stats::uniroot(function(t) eval_baseline(t, b_name, b_par)$H0 - u_val,
                         interval = c(1e-6, 1000), extendInt = "yes")$root,
          error = function(e) sqrt(u_val)
        )
        res
      })
    }
  }
  
  if (grepl("bivariate", type)) {
    u1 <- stats::runif(n)
    u2 <- stats::runif(n)
    
    if (type == "multiplicative_bivariate") {
      H1_val <- -log(u1) / (z * rho)
      H2_val <- -log(u2) / (z * rho)
      t1_true <- inv_H0(H1_val, baseline, baseline_params)
      t2_true <- inv_H0(H2_val, baseline, baseline_params)
    } else {
      t1_true <- sapply(1:n, function(i) {
        val <- -log(u1[i])
        stats::uniroot(function(t) eval_baseline(t, baseline, baseline_params)$H0 + z[i] * rho[i] * t - val,
                       interval = c(1e-6, 1000), extendInt = "yes")$root
      })
      t2_true <- sapply(1:n, function(i) {
        val <- -log(u2[i])
        stats::uniroot(function(t) eval_baseline(t, baseline, baseline_params)$H0 + z[i] * rho[i] * t - val,
                       interval = c(1e-6, 1000), extendInt = "yes")$root
      })
    }
    
    c_max1 <- stats::quantile(t1_true, 1 - censor_rate, names = FALSE)
    c_max2 <- stats::quantile(t2_true, 1 - censor_rate, names = FALSE)
    c1 <- stats::runif(n, 0, max(c_max1 * 1.5, 1))
    c2 <- stats::runif(n, 0, max(c_max2 * 1.5, 1))
    
    time1 <- pmin(t1_true, c1)
    status1 <- as.numeric(t1_true <= c1)
    time2 <- pmin(t2_true, c2)
    status2 <- as.numeric(t2_true <= c2)
    
    data.frame(id = 1:n, time1 = time1, time2 = time2,
               status1 = status1, status2 = status2,
               x1 = x1, x2 = x2, frailty = z)
  } else {
    u <- stats::runif(n)
    if (type == "multiplicative_univariate") {
      H_val <- -log(u) / (z * rho)
      t_true <- inv_H0(H_val, baseline, baseline_params)
    } else {
      t_true <- sapply(1:n, function(i) {
        val <- -log(u[i])
        stats::uniroot(function(t) eval_baseline(t, baseline, baseline_params)$H0 + z[i] * rho[i] * t - val,
                       interval = c(1e-6, 1000), extendInt = "yes")$root
      })
    }
    c_max <- stats::quantile(t_true, 1 - censor_rate, names = FALSE)
    c_time <- stats::runif(n, 0, max(c_max * 1.5, 1))
    time <- pmin(t_true, c_time)
    status <- as.numeric(t_true <= c_time)
    
    data.frame(id = 1:n, time = time, status = status,
               x1 = x1, x2 = x2, frailty = z)
  }
}
