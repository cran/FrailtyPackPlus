#' Correlated Frailty Log-Survival and Bivariate Derivatives
#'
#' @description
#' Core functions implementing bivariate Laplace transforms, log-survival functions,
#' and their partial derivatives with respect to cumulative baseline hazards s1 and s2.
#'
#' @name frailty_core
#' @keywords internal
NULL

frailty_param_names <- function(frailty) {
  frailty <- match.arg(frailty, c("gamma", "invgaussian", "posstable", "comppoisson", "none"))
  switch(frailty,
    gamma = c("sigma", "rho"),
    invgaussian = c("sigma", "rho"),
    posstable = c("theta", "rho"),
    comppoisson = c("eta", "nu", "theta"),
    none = character(0)
  )
}

frailty_param_count <- function(frailty) {
  length(frailty_param_names(frailty))
}

log_phi_derivatives <- function(s1, s2, frailty, params) {
  s1 <- pmax(as.numeric(s1), 1e-12)
  s2 <- pmax(as.numeric(s2), 1e-12)
  frailty <- match.arg(frailty, c("gamma", "invgaussian", "posstable", "comppoisson", "none"))
  
  n <- length(s1)
  
  if (frailty == "none") {
    phi <- -s1 - s2
    phi1 <- rep(-1.0, n)
    phi2 <- rep(-1.0, n)
    phi12 <- rep(0.0, n)
  } else if (frailty == "gamma") {
    sigma <- max(1e-6, params[1])
    rho <- min(max(1e-6, params[2]), 1 - 1e-6)
    sig2 <- sigma^2
    
    t1 <- 1 + sig2 * s1
    t2 <- 1 + sig2 * s2
    t12 <- 1 + sig2 * (s1 + s2)
    
    phi <- - (rho / sig2) * log(t12) - ((1 - rho) / sig2) * log(t1) - ((1 - rho) / sig2) * log(t2)
    phi1 <- - rho / t12 - (1 - rho) / t1
    phi2 <- - rho / t12 - (1 - rho) / t2
    phi12 <- rho * sig2 / (t12^2)
  } else if (frailty == "invgaussian") {
    sigma <- max(1e-6, params[1])
    rho <- min(max(1e-6, params[2]), 1 - 1e-6)
    sig2 <- sigma^2
    
    r1 <- sqrt(1 + 2 * sig2 * s1)
    r2 <- sqrt(1 + 2 * sig2 * s2)
    r12 <- sqrt(1 + 2 * sig2 * (s1 + s2))
    
    phi <- (rho / sig2) * (1 - r12) + ((1 - rho) / sig2) * (1 - r1) + ((1 - rho) / sig2) * (1 - r2)
    phi1 <- - rho / r12 - (1 - rho) / r1
    phi2 <- - rho / r12 - (1 - rho) / r2
    phi12 <- rho * sig2 / (r12^3)
  } else if (frailty == "posstable") {
    theta <- min(max(1e-4, params[1]), 1 - 1e-4)
    rho <- min(max(1e-6, params[2]), 1 - 1e-6)
    
    phi <- - rho * ((s1 + s2)^theta) - (1 - rho) * (s1^theta) - (1 - rho) * (s2^theta)
    phi1 <- - rho * theta * ((s1 + s2)^(theta - 1)) - (1 - rho) * theta * (s1^(theta - 1))
    phi2 <- - rho * theta * ((s1 + s2)^(theta - 1)) - (1 - rho) * theta * (s2^(theta - 1))
    phi12 <- rho * theta * (1 - theta) * ((s1 + s2)^(theta - 2))
  } else if (frailty == "comppoisson") {
    eta <- max(1e-4, params[1])
    nu <- max(1e-4, params[2])
    theta <- max(1e-4, params[3])
    
    t1 <- (nu / (nu + s1))^eta
    t2 <- (nu / (nu + s2))^eta
    A <- pmax(1e-12, theta + 2 - t1 - t2)
    
    phi <- (nu * theta / eta) * (log(theta) - log(A))
    phi1 <- - (nu * theta * (nu^eta)) / (A * ((nu + s1)^(eta + 1)))
    phi2 <- - (nu * theta * (nu^eta)) / (A * ((nu + s2)^(eta + 1)))
    phi12 <- (nu * theta * eta * (nu^(2 * eta))) / ((A^2) * ((nu + s1)^(eta + 1)) * ((nu + s2)^(eta + 1)))
  }
  
  # Ensure numerical validity
  phi <- pmin(phi, 0)
  phi1 <- pmin(phi1, -1e-12)
  phi2 <- pmin(phi2, -1e-12)
  phi12 <- pmax(phi12, 0)
  
  list(phi = phi, phi1 = phi1, phi2 = phi2, phi12 = phi12)
}

bivariate_loglik_contrib <- function(t1, d1, t2, d2, eta1, eta2,
                                     baseline1, base_params1,
                                     baseline2, base_params2,
                                     frailty, frailty_params) {
  res1 <- baseline_H_h(t1, baseline1, base_params1)
  res2 <- baseline_H_h(t2, baseline2, base_params2)
  
  s1 <- eta1 * res1$H
  s2 <- eta2 * res2$H
  
  pder <- log_phi_derivatives(s1, s2, frailty, frailty_params)
  
  log_h1 <- log(eta1) + log(res1$h)
  log_h2 <- log(eta2) + log(res2$h)
  
  ll <- numeric(length(t1))
  
  # Case 1: d1 == 1, d2 == 1 (both uncensored)
  idx11 <- which(d1 == 1 & d2 == 1)
  if (length(idx11) > 0) {
    cross_term <- pmax(1e-15, pder$phi12[idx11] + pder$phi1[idx11] * pder$phi2[idx11])
    ll[idx11] <- pder$phi[idx11] + log(cross_term) + log_h1[idx11] + log_h2[idx11]
  }
  
  # Case 2: d1 == 1, d2 == 0 (t1 uncensored, t2 censored)
  idx10 <- which(d1 == 1 & d2 == 0)
  if (length(idx10) > 0) {
    term1 <- pmax(1e-15, -pder$phi1[idx10])
    ll[idx10] <- pder$phi[idx10] + log(term1) + log_h1[idx10]
  }
  
  # Case 3: d1 == 0, d2 == 1 (t1 censored, t2 uncensored)
  idx01 <- which(d1 == 0 & d2 == 1)
  if (length(idx01) > 0) {
    term2 <- pmax(1e-15, -pder$phi2[idx01])
    ll[idx01] <- pder$phi[idx01] + log(term2) + log_h2[idx01]
  }
  
  # Case 4: d1 == 0, d2 == 0 (both censored)
  idx00 <- which(d1 == 0 & d2 == 0)
  if (length(idx00) > 0) {
    ll[idx00] <- pder$phi[idx00]
  }
  
  ll
}
