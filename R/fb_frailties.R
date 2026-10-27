#' Frailty Laplace Transform and Derivatives
#'
#' Evaluates Laplace transform L0(s) = E[exp(-s Z)], first derivative magnitude
#' L1(s) = -L'(s) = E[Z exp(-s Z)], and second derivative L2(s) = L''(s) = E[Z^2 exp(-s Z)]
#' for eight frailty distributions.
#'
#' @param s Non-negative numeric vector.
#' @param frailty Character specifying frailty distribution: \code{"gamma"},
#'   \code{"invgauss"}, \code{"posstable"}, \code{"compoundpoisson"},
#'   \code{"compoundnb"}, \code{"genlindley"}, \code{"weightedlindley"},
#'   or \code{"mixture_invgauss"}.
#' @param params Numeric vector of frailty parameters.
#'
#' @return A list with elements \code{L0}, \code{L1}, \code{L2}, \code{log_L0},
#'   \code{log_L1}, and \code{log_L2}.
#' @export
#' @keywords internal
eval_laplace <- function(s, frailty = c("gamma", "invgauss", "posstable",
                                        "compoundpoisson", "compoundnb",
                                        "genlindley", "weightedlindley",
                                        "mixture_invgauss"), params) {
  frailty <- match.arg(frailty)
  s <- pmax(as.numeric(s), 0)
  eps <- 1e-15
  
  if (frailty == "gamma") {
    theta <- max(params[1], 1e-8)
    u <- 1 + s / theta
    log_u <- log(u)
    log_L0 <- -theta * log_u
    log_L1 <- -(theta + 1) * log_u
    log_L2 <- log(theta + 1) - log(theta) - (theta + 2) * log_u
    L0 <- exp(log_L0)
    L1 <- exp(log_L1)
    L2 <- exp(log_L2)
  } else if (frailty == "invgauss") {
    theta <- max(params[1], 1e-8)
    w <- sqrt(1 + 2 * theta * s)
    log_L0 <- (1 - w) / theta
    log_L1 <- log_L0 - log(w)
    log_L2 <- log_L0 - 2 * log(w) + log(pmax(1 + theta / w, eps))
    L0 <- exp(log_L0)
    L1 <- exp(log_L1)
    L2 <- exp(log_L2)
  } else if (frailty == "posstable") {
    alpha <- min(max(params[1], 1e-6), 0.9999)
    s_pos <- pmax(s, 1e-12)
    s_a <- s_pos^alpha
    log_L0 <- -s_a
    log_L1 <- log(alpha) + (alpha - 1) * log(s_pos) - s_a
    term2 <- pmax(alpha * s_a + (1 - alpha), eps)
    log_L2 <- log(alpha) + (alpha - 2) * log(s_pos) + log(term2) - s_a
    L0 <- exp(log_L0)
    L1 <- exp(log_L1)
    L2 <- exp(log_L2)
  } else if (frailty == "compoundpoisson") {
    theta1 <- max(params[1], 1e-8) # Poisson mean parameter
    theta2 <- max(params[2], 1e-8) # shape parameter
    theta3 <- if (length(params) >= 3) max(params[3], 1e-8) else 1.0
    u <- 1 + theta3 * s
    log_u <- log(u)
    g <- (theta1 / theta2) * (1 - exp(-theta2 * log_u))
    log_L0 <- -g
    log_L1 <- log(theta1) + log(theta3) - (theta2 + 1) * log_u - g
    L0 <- exp(log_L0)
    L1 <- exp(log_L1)
    L2_term <- (theta1 * theta3 * exp(-(theta2 + 1) * log_u))^2 + 
               theta1 * (theta3^2) * (theta2 + 1) * exp(-(theta2 + 2) * log_u)
    L2 <- pmax(L2_term * L0, eps)
    log_L2 <- log(L2)
  } else if (frailty == "compoundnb") {
    theta1 <- max(params[1], 1e-8)
    theta2 <- max(params[2], 1e-8)
    theta3 <- if (length(params) >= 3) max(params[3], 1e-8) else 1.0
    alpha <- if (length(params) >= 4) max(params[4], 1e-8) else 1.0
    u <- 1 + theta3 * s
    log_u <- log(u)
    g <- (theta1 / theta2) * (1 - exp(-alpha * log_u))
    v <- pmax(1 + g, eps)
    log_v <- log(v)
    log_L0 <- -theta2 * log_v
    gp <- theta1 * theta3 * alpha * exp(-(alpha + 1) * log_u) / theta2
    gpp <- -theta1 * (theta3^2) * alpha * (alpha + 1) * exp(-(alpha + 2) * log_u) / theta2
    L0 <- exp(log_L0)
    L1 <- pmax(theta2 * exp(-(theta2 + 1) * log_v) * gp, eps)
    log_L1 <- log(L1)
    L2_term <- theta2 * (theta2 + 1) * exp(-(theta2 + 2) * log_v) * (gp^2) - 
               theta2 * exp(-(theta2 + 1) * log_v) * gpp
    L2 <- pmax(L2_term, eps)
    log_L2 <- log(L2)
  } else if (frailty == "genlindley") {
    theta <- max(params[1], 1e-8)
    alpha <- max(params[2], 1e-8)
    d <- theta + s
    log_C <- 2 * log(theta) - log(theta + alpha)
    log_L0 <- log_C + log(d + alpha) - 2 * log(d)
    log_L1 <- log_C + log(d + 2 * alpha) - 3 * log(d)
    log_L2 <- log(2) + log_C + log(d + 3 * alpha) - 4 * log(d)
    L0 <- exp(log_L0)
    L1 <- exp(log_L1)
    L2 <- exp(log_L2)
  } else if (frailty == "weightedlindley") {
    theta <- max(params[1], 1e-8)
    c_par <- max(params[2], 1e-8)
    d <- theta + s
    log_C <- c_par * log(theta) - log(c_par + theta)
    log_L0 <- log_C + log(c_par + theta + s) - c_par * log(d)
    log_L1 <- log_C + log(c_par * (c_par + theta + s) + d) - (c_par + 1) * log(d)
    term_l2 <- c_par * (c_par + 1) * (c_par + theta + s) + 2 * (c_par + 1) * d
    log_L2 <- log_C + log(pmax(term_l2, eps)) - (c_par + 2) * log(d)
    L0 <- exp(log_L0)
    L1 <- exp(log_L1)
    L2 <- exp(log_L2)
  } else if (frailty == "mixture_invgauss") {
    p1 <- min(max(params[1], 1e-6), 0.9999)
    theta1 <- max(params[2], 1e-8)
    theta2 <- max(params[3], 1e-8)
    p2 <- 1 - p1
    res1 <- eval_laplace(s, "invgauss", theta1)
    res2 <- eval_laplace(s, "invgauss", theta2)
    L0 <- p1 * res1$L0 + p2 * res2$L0
    L1 <- p1 * res1$L1 + p2 * res2$L1
    L2 <- p1 * res1$L2 + p2 * res2$L2
    log_L0 <- log(pmax(L0, eps))
    log_L1 <- log(pmax(L1, eps))
    log_L2 <- log(pmax(L2, eps))
  }
  
  list(L0 = L0, L1 = L1, L2 = L2, log_L0 = log_L0, log_L1 = log_L1, log_L2 = log_L2)
}

#' Random Generation for Frailty Distributions
#'
#' Simulates frailty random variates.
#'
#' @param n Integer number of observations.
#' @param frailty Character string specifying frailty distribution.
#' @param params Numeric parameter vector.
#'
#' @return A numeric vector of length \code{n}.
#' @export
#' @keywords internal
r_frailty <- function(n, frailty = c("gamma", "invgauss", "posstable",
                                     "compoundpoisson", "compoundnb",
                                     "genlindley", "weightedlindley",
                                     "mixture_invgauss"), params) {
  frailty <- match.arg(frailty)
  if (frailty == "gamma") {
    theta <- params[1]
    stats::rgamma(n, shape = theta, rate = theta)
  } else if (frailty == "invgauss") {
    theta <- params[1]
    # Standard inverse Gaussian with mean 1, shape lambda = 1/theta
    mu <- 1
    lambda_ig <- 1 / theta
    y <- stats::rnorm(n)^2
    x <- mu + (mu^2 * y) / (2 * lambda_ig) - (mu / (2 * lambda_ig)) * sqrt(4 * mu * lambda_ig * y + mu^2 * y^2)
    z <- stats::runif(n)
    ifelse(z <= mu / (mu + x), x, mu^2 / x)
  } else if (frailty == "posstable") {
    alpha <- params[1]
    # Chambers-Mallows-Stuck algorithm for positive stable
    u <- stats::runif(n, 0, pi)
    e <- stats::rexp(n)
    ((sin(alpha * u) / (sin(u)^alpha)) * (sin((1 - alpha) * u) / e)^(1 - alpha))^(1 / alpha)
  } else if (frailty == "compoundpoisson") {
    theta1 <- params[1]
    theta2 <- params[2]
    theta3 <- if (length(params) >= 3) params[3] else 1.0
    N <- stats::rpois(n, lambda = theta1)
    sapply(N, function(k) if (k == 0) 0.001 else sum(stats::rgamma(k, shape = theta2, scale = theta3)))
  } else if (frailty == "compoundnb") {
    theta1 <- params[1]
    theta2 <- params[2]
    theta3 <- if (length(params) >= 3) params[3] else 1.0
    alpha <- if (length(params) >= 4) params[4] else 1.0
    N <- stats::rnbinom(n, size = theta2, prob = 1 / (1 + theta1 / theta2))
    sapply(N, function(k) if (k == 0) 0.001 else sum(stats::rgamma(k, shape = alpha, scale = theta3)))
  } else if (frailty == "genlindley") {
    theta <- params[1]
    alpha <- params[2]
    p <- theta / (theta + alpha)
    u <- stats::runif(n)
    ifelse(u <= p, stats::rexp(n, rate = theta), stats::rgamma(n, shape = 2, rate = theta))
  } else if (frailty == "weightedlindley") {
    theta <- params[1]
    c_par <- params[2]
    p <- theta / (theta + c_par)
    u <- stats::runif(n)
    ifelse(u <= p, stats::rgamma(n, shape = c_par, rate = theta), stats::rgamma(n, shape = c_par + 1, rate = theta))
  } else if (frailty == "mixture_invgauss") {
    p1 <- params[1]
    theta1 <- params[2]
    theta2 <- params[3]
    u <- stats::runif(n)
    z1 <- r_frailty(n, "invgauss", theta1)
    z2 <- r_frailty(n, "invgauss", theta2)
    ifelse(u <= p1, z1, z2)
  }
}
