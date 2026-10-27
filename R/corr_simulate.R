#' Simulate Bivariate Survival Data with Correlated Frailty
#'
#' @description
#' Generates synthetic bivariate right-censored survival data from correlated
#' frailty models with user-specified frailty families, baseline distributions,
#' regression coefficients, and censoring rates.
#'
#' @param n Integer, sample size (number of pairs/clusters). Default is 100.
#' @param frailty Character, frailty distribution. One of \code{"gamma"},
#'   \code{"invgaussian"}, \code{"posstable"}, \code{"comppoisson"}, or \code{"none"}.
#' @param frailty.params List of frailty parameters. E.g., for \code{"gamma"} or
#'   \code{"invgaussian"}: \code{list(sigma = 0.5, rho = 0.5)}; for \code{"posstable"}:
#'   \code{list(theta = 0.7, rho = 0.5)}; for \code{"comppoisson"}:
#'   \code{list(eta = 0.5, nu = 1.0, theta = 1.0)}.
#' @param baseline Character, baseline distribution for T1 and T2. One of
#'   \code{"genweibull"}, \code{"genloglogistic"}, \code{"weibull"},
#'   \code{"exponential"}, \code{"pareto"}, \code{"lfr"}, or \code{"burr"}.
#' @param baseline.params List of baseline parameters. E.g., for \code{"weibull"}:
#'   \code{list(lambda = 0.05, gamma = 1.2)}; for \code{"genweibull"} or
#'   \code{"genloglogistic"}: \code{list(alpha = 1.5, lambda = 0.05, gamma = 1.2)}.
#' @param beta Numeric vector of regression coefficients. Default is \code{c(0.5, -0.3)}.
#' @param censor.rate Numeric between 0 and 1, target censoring proportion. Default is 0.2.
#' @param seed Optional integer for reproducibility.
#'
#' @return A \code{data.frame} containing:
#' \item{id}{Integer cluster / pair identifier (1 to n).}
#' \item{t1}{Observed time for first component.}
#' \item{d1}{Censoring indicator for first component (1 = event, 0 = censored).}
#' \item{t2}{Observed time for second component.}
#' \item{d2}{Censoring indicator for second component (1 = event, 0 = censored).}
#' \item{x1, x2, ...}{Simulated covariates.}
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
#' @author
#' Shikhar Tyagi, Arvind Pandey, Vrijesh Tripathi
#'
#' @export
#' @examples
#' dat <- sim_corrFrailty(n = 50, frailty = "gamma",
#'                        frailty.params = list(sigma = 0.5, rho = 0.4),
#'                        baseline = "weibull",
#'                        baseline.params = list(lambda = 0.05, gamma = 1.2),
#'                        beta = c(0.4, -0.3), seed = 123)
#' head(dat)
sim_corrFrailty <- function(n = 100,
                            frailty = c("gamma", "invgaussian", "posstable", "comppoisson", "none"),
                            frailty.params = list(sigma = 0.5, rho = 0.5),
                            baseline = c("genweibull", "genloglogistic", "weibull",
                                         "exponential", "pareto", "lfr", "burr"),
                            baseline.params = list(lambda = 0.05, gamma = 1.2),
                            beta = c(0.5, -0.3),
                            censor.rate = 0.2,
                            seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  frailty <- match.arg(frailty)
  baseline <- match.arg(baseline)
  
  # Simulate covariates: X1 ~ Bernoulli(0.5), X2 ~ Normal(0, 1)
  p_cov <- length(beta)
  X <- if (p_cov > 0) {
    x1 <- stats::rbinom(n, 1, 0.5)
    if (p_cov == 1) {
      matrix(x1, ncol = 1, dimnames = list(NULL, "x1"))
    } else {
      x_rest <- matrix(stats::rnorm(n * (p_cov - 1)), nrow = n, ncol = p_cov - 1)
      colnames(x_rest) <- paste0("x", 2:p_cov)
      cbind(x1 = x1, x_rest)
    }
  } else {
    NULL
  }
  
  eta <- if (p_cov > 0) as.vector(exp(X %*% beta)) else rep(1, n)
  
  # Simulate Frailties Z1, Z2
  if (frailty == "none") {
    Z1 <- rep(1, n)
    Z2 <- rep(1, n)
  } else if (frailty == "gamma") {
    sigma <- frailty.params$sigma
    rho <- frailty.params$rho
    sig2 <- sigma^2
    
    # Bivariate Gamma representation: Z1 = V0 + V1, Z2 = V0 + V2
    # V0 ~ Gamma(rho/sig2, 1/sig2), V1, V2 ~ Gamma((1-rho)/sig2, 1/sig2)
    k0 <- rho / sig2
    k1 <- (1 - rho) / sig2
    scale_val <- sig2
    
    V0 <- stats::rgamma(n, shape = max(1e-4, k0), scale = scale_val)
    V1 <- stats::rgamma(n, shape = max(1e-4, k1), scale = scale_val)
    V2 <- stats::rgamma(n, shape = max(1e-4, k1), scale = scale_val)
    
    Z1 <- V0 + V1
    Z2 <- V0 + V2
  } else if (frailty == "invgaussian") {
    sigma <- frailty.params$sigma
    rho <- frailty.params$rho
    sig2 <- sigma^2
    
    # Simulate correlated IG: common shock + independent IG shocks
    # Standard IG sampler
    rinvgauss <- function(n_samp, mu, lambda_param) {
      y <- stats::rnorm(n_samp)^2
      x <- mu + (mu^2 * y) / (2 * lambda_param) -
        (mu / (2 * lambda_param)) * sqrt(4 * mu * lambda_param * y + mu^2 * y^2)
      u <- stats::runif(n_samp)
      ifelse(u <= mu / (mu + x), x, mu^2 / x)
    }
    
    lambda_0 <- rho / sig2
    lambda_1 <- (1 - rho) / sig2
    
    # Common and individual shocks
    W0 <- rinvgauss(n, mu = rho, lambda_param = max(1e-4, lambda_0))
    W1 <- rinvgauss(n, mu = 1 - rho, lambda_param = max(1e-4, lambda_1))
    W2 <- rinvgauss(n, mu = 1 - rho, lambda_param = max(1e-4, lambda_1))
    
    Z1 <- W0 + W1
    Z2 <- W0 + W2
  } else if (frailty == "posstable") {
    theta <- frailty.params$theta
    rho <- frailty.params$rho
    
    # Positive stable generation via Chambers-Mallows-Stuck algorithm
    rposstable <- function(n_samp, alpha_val) {
      u <- stats::runif(n_samp, -pi / 2, pi / 2)
      w <- stats::rexp(n_samp)
      (sin(alpha_val * (u + pi / 2)) / (cos(u))^(1 / alpha_val)) *
        ((cos(u - alpha_val * (u + pi / 2)) / w)^((1 - alpha_val) / alpha_val))
    }
    
    S0 <- rposstable(n, theta) * (rho^(1 / theta))
    S1 <- rposstable(n, theta) * ((1 - rho)^(1 / theta))
    S2 <- rposstable(n, theta) * ((1 - rho)^(1 / theta))
    
    Z1 <- S0 + S1
    Z2 <- S0 + S2
  } else if (frailty == "comppoisson") {
    eta_param <- frailty.params$eta
    nu_param <- frailty.params$nu
    theta_param <- frailty.params$theta
    
    # Compound Poisson: Poisson(lambda_CP) sum of Gamma
    lambda_CP <- stats::rgamma(n, shape = theta_param, rate = 1)
    N1 <- stats::rpois(n, lambda_CP)
    N2 <- stats::rpois(n, lambda_CP)
    
    Z1 <- vapply(N1, function(k) if (k == 0) 0.01 else sum(stats::rgamma(k, shape = eta_param, rate = nu_param)), numeric(1))
    Z2 <- vapply(N2, function(k) if (k == 0) 0.01 else sum(stats::rgamma(k, shape = eta_param, rate = nu_param)), numeric(1))
  }
  
  Z1 <- pmax(1e-6, Z1)
  Z2 <- pmax(1e-6, Z2)
  
  # Generate Conditional True Lifetimes T1, T2:
  # S(t | Z, X) = exp(-Z * eta * H_0(t)) = U ~ Uniform(0, 1)
  # -> H_0(T) = -log(U) / (Z * eta) -> T = H_0^{-1}(-log(U) / (Z * eta))
  U1 <- stats::runif(n)
  U2 <- stats::runif(n)
  
  req_H1 <- -log(U1) / (Z1 * eta)
  req_H2 <- -log(U2) / (Z2 * eta)
  
  b_vec <- unlist(baseline.params[baseline_param_names(baseline)])
  T1_true <- baseline_inv_H(req_H1, baseline, b_vec)
  T2_true <- baseline_inv_H(req_H2, baseline, b_vec)
  
  # Generate Censoring Times
  if (censor.rate > 0) {
    # Exponential censoring to achieve target rate
    rate_c1 <- censor.rate / (1 - censor.rate + 1e-6) / mean(T1_true)
    rate_c2 <- censor.rate / (1 - censor.rate + 1e-6) / mean(T2_true)
    
    C1 <- stats::rexp(n, rate = max(1e-6, rate_c1))
    C2 <- stats::rexp(n, rate = max(1e-6, rate_c2))
    
    t1_obs <- pmin(T1_true, C1)
    d1_obs <- as.integer(T1_true <= C1)
    
    t2_obs <- pmin(T2_true, C2)
    d2_obs <- as.integer(T2_true <= C2)
  } else {
    t1_obs <- T1_true
    d1_obs <- rep(1L, n)
    t2_obs <- T2_true
    d2_obs <- rep(1L, n)
  }
  
  res_df <- data.frame(
    id = seq_len(n),
    t1 = round(t1_obs, 4),
    d1 = d1_obs,
    t2 = round(t2_obs, 4),
    d2 = d2_obs
  )
  
  if (p_cov > 0) {
    res_df <- cbind(res_df, X)
  }
  
  res_df
}
