# Distribution Functions, Laplace Transforms, and Dependence Measures
# Shared Frailty Models under Reversed Hazard Rate (RHR)

#' Baseline Probability Density Function under Reversed Hazard Rate
#'
#' Computes the probability density function (PDF) for various baseline distributions
#' used in reversed hazard rate frailty models.
#'
#' @param t Numeric vector of positive failure times.
#' @param dist Character string specifying baseline distribution. One of \code{"GE"},
#'   \code{"EG"}, \code{"GIR"}, \code{"MIW"}, \code{"GLL-I"}, \code{"GLL-II"},
#'   \code{"GIE"}, \code{"Weibull"}, \code{"Exponential"}, \code{"Log-Logistic"}.
#' @param params Numeric named vector of baseline parameters.
#'
#' @return Numeric vector of density values.
#' @export
#'
#' @examples
#' d_rhr_baseline(c(10, 20, 30), dist = "GE", params = c(zeta = 12.4, delta = 0.11))
#' d_rhr_baseline(c(10, 20, 30), dist = "EG", params = c(zeta = 9.1, delta = 0.10))
d_rhr_baseline <- function(t, dist = c("GE", "EG", "GIR", "MIW", "GLL-I", "GLL-II", "GIE", "Weibull", "Exponential", "Log-Logistic"), params) {
  dist <- match.arg(dist)
  t <- pmax(t, 1e-10)
  
  if (dist == "GE") {
    zeta <- as.numeric(params["zeta"])
    delta <- as.numeric(params["delta"])
    ed <- exp(-delta * t)
    f <- zeta * delta * ed * (pmax(1 - ed, 1e-12))^(zeta - 1)
  } else if (dist == "EG") {
    zeta <- as.numeric(params["zeta"])
    delta <- as.numeric(params["delta"])
    ed <- exp(-delta * t)
    f <- zeta * delta * ed * exp(-zeta * ed)
  } else if (dist == "GIR") {
    zeta <- as.numeric(params["zeta"])
    delta <- as.numeric(params["delta"])
    z <- exp(-delta / (t^2))
    f <- (2 * delta * zeta / (t^3)) * z * (pmax(1 - z, 1e-12))^(zeta - 1)
  } else if (dist == "MIW") {
    zeta <- as.numeric(params["zeta"])
    delta <- as.numeric(params["delta"])
    xi <- as.numeric(params["xi"])
    term <- zeta * (t^(-delta)) * exp(-xi * t)
    psi0 <- term * (delta / t + xi)
    f <- psi0 * exp(-term)
  } else if (dist == "GLL-I") {
    lambda <- as.numeric(params["lambda"])
    alpha <- as.numeric(params["alpha"])
    gamma_val <- as.numeric(params["gamma"])
    ltg <- (lambda * t)^gamma_val
    F_val <- (ltg / (1 + ltg))^alpha
    r_val <- (alpha * gamma_val) / (t * (1 + ltg))
    f <- r_val * F_val
  } else if (dist == "GLL-II") {
    lambda <- as.numeric(params["lambda"])
    alpha <- as.numeric(params["alpha"])
    gamma_val <- as.numeric(params["gamma"])
    ltg <- (lambda * t)^gamma_val
    f <- alpha * gamma_val * (lambda^gamma_val) * (t^(gamma_val - 1)) * (1 + ltg)^(-alpha - 1)
  } else if (dist == "GIE") {
    zeta <- as.numeric(params["zeta"])
    delta <- as.numeric(params["delta"])
    z <- exp(-delta / t)
    f <- (zeta * delta / (t^2)) * z * (pmax(1 - z, 1e-12))^(zeta - 1)
  } else if (dist == "Weibull") {
    zeta <- as.numeric(params["zeta"])
    delta <- as.numeric(params["delta"])
    f <- stats::dweibull(t, shape = zeta, scale = 1 / delta)
  } else if (dist == "Exponential") {
    delta <- as.numeric(params["delta"])
    f <- stats::dexp(t, rate = delta)
  } else if (dist == "Log-Logistic") {
    lambda <- as.numeric(params["lambda"])
    gamma_val <- as.numeric(params["gamma"])
    ltg <- (lambda * t)^gamma_val
    f <- (gamma_val * (lambda^gamma_val) * (t^(gamma_val - 1))) / ((1 + ltg)^2)
  } else {
    stop("Unknown baseline distribution: ", dist)
  }
  return(pmax(f, 0))
}

#' Baseline Cumulative Distribution Function under Reversed Hazard Rate
#'
#' Computes the cumulative distribution function (CDF) for baseline distributions
#' under reversed hazard rate models.
#'
#' @param t Numeric vector of positive failure times.
#' @param dist Character string specifying baseline distribution. One of \code{"GE"},
#'   \code{"EG"}, \code{"GIR"}, \code{"MIW"}, \code{"GLL-I"}, \code{"GLL-II"},
#'   \code{"GIE"}, \code{"Weibull"}, \code{"Exponential"}, \code{"Log-Logistic"}.
#' @param params Numeric named vector of baseline parameters.
#'
#' @return Numeric vector of cumulative probabilities in (0, 1).
#' @export
#'
#' @examples
#' p_rhr_baseline(c(10, 20, 30), dist = "GE", params = c(zeta = 12.4, delta = 0.11))
p_rhr_baseline <- function(t, dist = c("GE", "EG", "GIR", "MIW", "GLL-I", "GLL-II", "GIE", "Weibull", "Exponential", "Log-Logistic"), params) {
  dist <- match.arg(dist)
  t <- pmax(t, 1e-10)
  
  if (dist == "GE") {
    zeta <- as.numeric(params["zeta"])
    delta <- as.numeric(params["delta"])
    ed <- exp(-delta * t)
    F_val <- (pmax(1 - ed, 1e-12))^zeta
  } else if (dist == "EG") {
    zeta <- as.numeric(params["zeta"])
    delta <- as.numeric(params["delta"])
    ed <- exp(-delta * t)
    F_val <- exp(-zeta * ed)
  } else if (dist == "GIR") {
    zeta <- as.numeric(params["zeta"])
    delta <- as.numeric(params["delta"])
    z <- exp(-delta / (t^2))
    F_val <- 1 - (pmax(1 - z, 1e-12))^zeta
  } else if (dist == "MIW") {
    zeta <- as.numeric(params["zeta"])
    delta <- as.numeric(params["delta"])
    xi <- as.numeric(params["xi"])
    F_val <- exp(-zeta * (t^(-delta)) * exp(-xi * t))
  } else if (dist == "GLL-I") {
    lambda <- as.numeric(params["lambda"])
    alpha <- as.numeric(params["alpha"])
    gamma_val <- as.numeric(params["gamma"])
    ltg <- (lambda * t)^gamma_val
    F_val <- (ltg / (1 + ltg))^alpha
  } else if (dist == "GLL-II") {
    lambda <- as.numeric(params["lambda"])
    alpha <- as.numeric(params["alpha"])
    gamma_val <- as.numeric(params["gamma"])
    ltg <- (lambda * t)^gamma_val
    F_val <- 1 - (1 + ltg)^(-alpha)
  } else if (dist == "GIE") {
    zeta <- as.numeric(params["zeta"])
    delta <- as.numeric(params["delta"])
    z <- exp(-delta / t)
    F_val <- 1 - (pmax(1 - z, 1e-12))^zeta
  } else if (dist == "Weibull") {
    zeta <- as.numeric(params["zeta"])
    delta <- as.numeric(params["delta"])
    F_val <- stats::pweibull(t, shape = zeta, scale = 1 / delta)
  } else if (dist == "Exponential") {
    delta <- as.numeric(params["delta"])
    F_val <- stats::pexp(t, rate = delta)
  } else if (dist == "Log-Logistic") {
    lambda <- as.numeric(params["lambda"])
    gamma_val <- as.numeric(params["gamma"])
    ltg <- (lambda * t)^gamma_val
    F_val <- ltg / (1 + ltg)
  } else {
    stop("Unknown baseline distribution: ", dist)
  }
  return(pmin(pmax(F_val, 1e-15), 1 - 1e-15))
}

#' Baseline Reversed Hazard Rate Function
#'
#' Computes the reversed hazard rate (RHR), defined as \eqn{r_0(t) = f_0(t) / F_0(t)}.
#'
#' @param t Numeric vector of positive failure times.
#' @param dist Character string specifying baseline distribution.
#' @param params Numeric named vector of baseline parameters.
#'
#' @return Numeric vector of reversed hazard rate values.
#' @export
#'
#' @examples
#' rhr_baseline(c(10, 20, 30), dist = "GE", params = c(zeta = 12.4, delta = 0.11))
rhr_baseline <- function(t, dist = c("GE", "EG", "GIR", "MIW", "GLL-I", "GLL-II", "GIE", "Weibull", "Exponential", "Log-Logistic"), params) {
  dist <- match.arg(dist)
  t <- pmax(t, 1e-10)
  
  if (dist == "GE") {
    zeta <- as.numeric(params["zeta"])
    delta <- as.numeric(params["delta"])
    ed <- exp(-delta * t)
    r <- (zeta * delta * ed) / (pmax(1 - ed, 1e-12))
  } else if (dist == "EG") {
    zeta <- as.numeric(params["zeta"])
    delta <- as.numeric(params["delta"])
    r <- zeta * delta * exp(-delta * t)
  } else if (dist == "GIR") {
    zeta <- as.numeric(params["zeta"])
    delta <- as.numeric(params["delta"])
    z <- exp(-delta / (t^2))
    f <- (2 * delta * zeta / (t^3)) * z * (pmax(1 - z, 1e-12))^(zeta - 1)
    F_val <- pmax(1 - (pmax(1 - z, 1e-12))^zeta, 1e-12)
    r <- f / F_val
  } else if (dist == "MIW") {
    zeta <- as.numeric(params["zeta"])
    delta <- as.numeric(params["delta"])
    xi <- as.numeric(params["xi"])
    r <- zeta * (t^(-1 - delta)) * exp(-xi * t) * (delta + xi * t)
  } else if (dist == "GLL-I") {
    alpha <- as.numeric(params["alpha"])
    lambda <- as.numeric(params["lambda"])
    gamma_val <- as.numeric(params["gamma"])
    ltg <- (lambda * t)^gamma_val
    r <- (alpha * gamma_val) / (t * (1 + ltg))
  } else if (dist == "GLL-II") {
    alpha <- as.numeric(params["alpha"])
    lambda <- as.numeric(params["lambda"])
    gamma_val <- as.numeric(params["gamma"])
    ltg <- (lambda * t)^gamma_val
    f <- alpha * gamma_val * (lambda^gamma_val) * (t^(gamma_val - 1)) * (1 + ltg)^(-alpha - 1)
    F_val <- pmax(1 - (1 + ltg)^(-alpha), 1e-12)
    r <- f / F_val
  } else if (dist == "GIE") {
    zeta <- as.numeric(params["zeta"])
    delta <- as.numeric(params["delta"])
    z <- exp(-delta / t)
    f <- (zeta * delta / (t^2)) * z * (pmax(1 - z, 1e-12))^(zeta - 1)
    F_val <- pmax(1 - (pmax(1 - z, 1e-12))^zeta, 1e-12)
    r <- f / F_val
  } else {
    f <- d_rhr_baseline(t, dist, params)
    F_val <- p_rhr_baseline(t, dist, params)
    r <- f / pmax(F_val, 1e-12)
  }
  return(pmax(r, 1e-12))
}

#' Baseline Cumulative Reversed Hazard Rate Function
#'
#' Computes the cumulative reversed hazard rate (CRHR), defined as \eqn{R_0(t) = -\ln F_0(t)}.
#'
#' @param t Numeric vector of positive failure times.
#' @param dist Character string specifying baseline distribution.
#' @param params Numeric named vector of baseline parameters.
#'
#' @return Numeric vector of cumulative reversed hazard rate values.
#' @export
#'
#' @examples
#' crhr_baseline(c(10, 20, 30), dist = "GE", params = c(zeta = 12.4, delta = 0.11))
crhr_baseline <- function(t, dist = c("GE", "EG", "GIR", "MIW", "GLL-I", "GLL-II", "GIE", "Weibull", "Exponential", "Log-Logistic"), params) {
  dist <- match.arg(dist)
  t <- pmax(t, 1e-10)
  
  if (dist == "GE") {
    zeta <- as.numeric(params["zeta"])
    delta <- as.numeric(params["delta"])
    ed <- exp(-delta * t)
    R_val <- -zeta * log(pmax(1 - ed, 1e-12))
  } else if (dist == "EG") {
    zeta <- as.numeric(params["zeta"])
    delta <- as.numeric(params["delta"])
    R_val <- zeta * exp(-delta * t)
  } else if (dist == "GIR") {
    zeta <- as.numeric(params["zeta"])
    delta <- as.numeric(params["delta"])
    z <- exp(-delta / (t^2))
    R_val <- -log(pmax(1 - (pmax(1 - z, 1e-12))^zeta, 1e-12))
  } else if (dist == "MIW") {
    zeta <- as.numeric(params["zeta"])
    delta <- as.numeric(params["delta"])
    xi <- as.numeric(params["xi"])
    R_val <- zeta * (t^(-delta)) * exp(-xi * t)
  } else if (dist == "GLL-I") {
    alpha <- as.numeric(params["alpha"])
    lambda <- as.numeric(params["lambda"])
    gamma_val <- as.numeric(params["gamma"])
    ltg <- (lambda * t)^gamma_val
    R_val <- -alpha * log(ltg / (1 + ltg))
  } else if (dist == "GLL-II") {
    alpha <- as.numeric(params["alpha"])
    lambda <- as.numeric(params["lambda"])
    gamma_val <- as.numeric(params["gamma"])
    ltg <- (lambda * t)^gamma_val
    R_val <- -log(pmax(1 - (1 + ltg)^(-alpha), 1e-12))
  } else if (dist == "GIE") {
    zeta <- as.numeric(params["zeta"])
    delta <- as.numeric(params["delta"])
    z <- exp(-delta / t)
    R_val <- -log(pmax(1 - (pmax(1 - z, 1e-12))^zeta, 1e-12))
  } else {
    F_val <- p_rhr_baseline(t, dist, params)
    R_val <- -log(pmax(F_val, 1e-12))
  }
  return(pmax(R_val, 1e-12))
}

#' Baseline Quantile Function under Reversed Hazard Rate
#'
#' Computes quantiles for baseline distributions used in reversed hazard rate models.
#'
#' @param p Numeric vector of probabilities in (0, 1).
#' @param dist Character string specifying baseline distribution.
#' @param params Numeric named vector of baseline parameters.
#'
#' @return Numeric vector of quantiles.
#' @export
#'
#' @examples
#' q_rhr_baseline(c(0.1, 0.5, 0.9), dist = "GE", params = c(zeta = 12.4, delta = 0.11))
q_rhr_baseline <- function(p, dist = c("GE", "EG", "GIR", "MIW", "GLL-I", "GLL-II", "GIE", "Weibull", "Exponential", "Log-Logistic"), params) {
  dist <- match.arg(dist)
  p <- pmin(pmax(p, 1e-8), 1 - 1e-8)
  
  if (dist == "GE") {
    zeta <- as.numeric(params["zeta"])
    delta <- as.numeric(params["delta"])
    q <- -1 / delta * log(pmax(1e-12, 1 - p^(1 / zeta)))
  } else if (dist == "EG") {
    zeta <- as.numeric(params["zeta"])
    delta <- as.numeric(params["delta"])
    q <- -1 / delta * log(pmax(1e-12, -log(p) / zeta))
  } else if (dist == "GIR") {
    zeta <- as.numeric(params["zeta"])
    delta <- as.numeric(params["delta"])
    val <- -log(pmax(1e-12, 1 - (1 - p)^(1 / zeta)))
    q <- sqrt(delta / pmax(val, 1e-12))
  } else if (dist == "GIE") {
    zeta <- as.numeric(params["zeta"])
    delta <- as.numeric(params["delta"])
    val <- -log(pmax(1e-12, 1 - (1 - p)^(1 / zeta)))
    q <- delta / pmax(val, 1e-12)
  } else if (dist == "GLL-I") {
    alpha <- as.numeric(params["alpha"])
    lambda <- as.numeric(params["lambda"])
    gamma_val <- as.numeric(params["gamma"])
    p_a <- p^(1 / alpha)
    q <- (1 / lambda) * (p_a / pmax(1 - p_a, 1e-12))^(1 / gamma_val)
  } else if (dist == "GLL-II") {
    alpha <- as.numeric(params["alpha"])
    lambda <- as.numeric(params["lambda"])
    gamma_val <- as.numeric(params["gamma"])
    term <- (1 - p)^(-1 / alpha) - 1
    q <- (1 / lambda) * (pmax(term, 1e-12))^(1 / gamma_val)
  } else if (dist == "Weibull") {
    zeta <- as.numeric(params["zeta"])
    delta <- as.numeric(params["delta"])
    q <- stats::qweibull(p, shape = zeta, scale = 1 / delta)
  } else if (dist == "Exponential") {
    delta <- as.numeric(params["delta"])
    q <- stats::qexp(p, rate = delta)
  } else if (dist == "Log-Logistic") {
    lambda <- as.numeric(params["lambda"])
    gamma_val <- as.numeric(params["gamma"])
    q <- (1 / lambda) * (p / pmax(1 - p, 1e-12))^(1 / gamma_val)
  } else {
    q <- vapply(p, function(prob) {
      obj <- function(x) p_rhr_baseline(x, dist, params) - prob
      res <- try(stats::uniroot(obj, interval = c(1e-4, 500), extendInt = "yes"), silent = TRUE)
      if (inherits(res, "try-error")) return(20.0)
      res$root
    }, numeric(1))
  }
  return(pmax(q, 1e-6))
}

#' Laplace Transform of Frailty Distributions
#'
#' Computes the Laplace transform \eqn{L_W(s) = E[e^{-sW}]} for various frailty distributions.
#'
#' @param s Numeric vector of non-negative arguments.
#' @param frailty Character string indicating the frailty distribution. One of \code{"WL"},
#'   \code{"GL"}, \code{"GIG"}, \code{"Gamma"}, \code{"IG"}, \code{"Lindley"}, \code{"None"}.
#' @param theta Numeric parameter for frailty distribution.
#' @param mu Optional secondary parameter for \code{"GL"}.
#' @param theta2 Optional second parameter for \code{"GIG"}.
#'
#' @return Numeric vector of Laplace transform values.
#' @export
#'
#' @examples
#' laplace_frailty(c(0.1, 0.5, 1.0), frailty = "WL", theta = 3.4)
laplace_frailty <- function(s, frailty = c("WL", "GL", "GIG", "Gamma", "IG", "Lindley", "None"), theta, mu = NULL, theta2 = NULL) {
  frailty <- match.arg(frailty)
  s <- pmax(s, 0)
  
  if (frailty == "WL") {
    C <- sqrt(theta * (theta + 1))
    L <- ((C + s + theta) / (C + theta)) * (C / (s + C))^(theta + 1)
  } else if (frailty == "GL") {
    if (is.null(mu)) mu <- 0.5 * (1 + theta)
    eta <- theta * (1 + theta - mu)
    L <- (1 / (1 + theta)) * ( (theta^(mu + 1)) / ((s + theta)^mu) + (theta^eta) / ((s + theta)^eta) )
  } else if (frailty == "GIG") {
    th2 <- if (is.null(theta2)) theta else theta2
    w1 <- theta / (theta + th2)
    w2 <- th2 / (theta + th2)
    L1 <- exp(theta * (1 - sqrt(1 + 2 * s / theta)))
    L2 <- exp(th2 * (1 - sqrt(1 + 2 * s / th2)))
    L <- w1 * L1 + w2 * L2
  } else if (frailty == "Gamma") {
    L <- (1 + s / theta)^(-theta)
  } else if (frailty == "IG") {
    L <- exp(theta * (1 - sqrt(1 + 2 * s / theta)))
  } else if (frailty == "Lindley") {
    L <- (theta^2 * (1 + theta + s)) / ((1 + theta) * (s + theta)^2)
  } else if (frailty == "None") {
    L <- exp(-s)
  }
  return(pmin(pmax(L, 1e-15), 1.0))
}

#' First Derivative of Laplace Transform of Frailty Distributions
#'
#' Computes \eqn{L_W'(s) = d L_W(s) / ds}. Note that \eqn{L_W'(s) \le 0}.
#'
#' @param s Numeric vector of non-negative arguments.
#' @param frailty Character string indicating the frailty distribution.
#' @param theta Numeric parameter for frailty distribution.
#' @param mu Optional secondary parameter for \code{"GL"}.
#' @param theta2 Optional second parameter for \code{"GIG"}.
#'
#' @return Numeric vector of first derivative values.
#' @export
#'
#' @examples
#' dlaplace_frailty(c(0.1, 0.5, 1.0), frailty = "WL", theta = 3.4)
dlaplace_frailty <- function(s, frailty = c("WL", "GL", "GIG", "Gamma", "IG", "Lindley", "None"), theta, mu = NULL, theta2 = NULL) {
  frailty <- match.arg(frailty)
  s <- pmax(s, 0)
  
  if (frailty == "WL") {
    C <- sqrt(theta * (theta + 1))
    g <- (C / (s + C))^(theta + 1)
    num_term <- -(theta * (s + C) + (theta + 1) * theta)
    den_term <- (C + theta) * (s + C)
    dL <- g * num_term / den_term
  } else if (frailty == "GL") {
    if (is.null(mu)) mu <- 0.5 * (1 + theta)
    eta <- theta * (1 + theta - mu)
    dL <- (1 / (1 + theta)) * ( -mu * (theta^(mu + 1)) / ((s + theta)^(mu + 1)) - eta * (theta^eta) / ((s + theta)^(eta + 1)) )
  } else if (frailty == "GIG") {
    th2 <- if (is.null(theta2)) theta else theta2
    w1 <- theta / (theta + th2)
    w2 <- th2 / (theta + th2)
    dL1 <- -(1 + 2 * s / theta)^(-0.5) * exp(theta * (1 - sqrt(1 + 2 * s / theta)))
    dL2 <- -(1 + 2 * s / th2)^(-0.5) * exp(th2 * (1 - sqrt(1 + 2 * s / th2)))
    dL <- w1 * dL1 + w2 * dL2
  } else if (frailty == "Gamma") {
    dL <- -(1 + s / theta)^(-theta - 1)
  } else if (frailty == "IG") {
    dL <- -(1 + 2 * s / theta)^(-0.5) * exp(theta * (1 - sqrt(1 + 2 * s / theta)))
  } else if (frailty == "Lindley") {
    dL <- - (theta^2 * (s + theta + 2)) / ((1 + theta) * (s + theta)^3)
  } else if (frailty == "None") {
    dL <- -exp(-s)
  }
  return(pmin(dL, -1e-15))
}

#' Second Derivative of Laplace Transform of Frailty Distributions
#'
#' Computes \eqn{L_W''(s) = d^2 L_W(s) / ds^2}. Note that \eqn{L_W''(s) > 0}.
#'
#' @param s Numeric vector of non-negative arguments.
#' @param frailty Character string indicating the frailty distribution.
#' @param theta Numeric parameter for frailty distribution.
#' @param mu Optional secondary parameter for \code{"GL"}.
#' @param theta2 Optional second parameter for \code{"GIG"}.
#'
#' @return Numeric vector of second derivative values.
#' @export
#'
#' @examples
#' d2laplace_frailty(c(0.1, 0.5, 1.0), frailty = "WL", theta = 3.4)
d2laplace_frailty <- function(s, frailty = c("WL", "GL", "GIG", "Gamma", "IG", "Lindley", "None"), theta, mu = NULL, theta2 = NULL) {
  frailty <- match.arg(frailty)
  s <- pmax(s, 0)
  
  if (frailty == "WL") {
    C <- sqrt(theta * (theta + 1))
    g <- (C / (s + C))^(theta + 1)
    num_term <- -(theta * (s + C) + (theta + 1) * theta)
    den_term <- (C + theta) * (s + C)
    term1 <- -(theta + 1) / (s + C) * (num_term / den_term)
    term2 <- (-theta * den_term - num_term * (C + theta)) / (den_term^2)
    d2L <- g * (term1 + term2)
  } else if (frailty == "GL") {
    if (is.null(mu)) mu <- 0.5 * (1 + theta)
    eta <- theta * (1 + theta - mu)
    d2L <- (1 / (1 + theta)) * ( mu * (mu + 1) * (theta^(mu + 1)) / ((s + theta)^(mu + 2)) + eta * (eta + 1) * (theta^eta) / ((s + theta)^(eta + 2)) )
  } else if (frailty == "GIG") {
    th2 <- if (is.null(theta2)) theta else theta2
    w1 <- theta / (theta + th2)
    w2 <- th2 / (theta + th2)
    d2L1 <- ((1 + 2 * s / theta)^(-1.0) + (1 / theta) * (1 + 2 * s / theta)^(-1.5)) * exp(theta * (1 - sqrt(1 + 2 * s / theta)))
    d2L2 <- ((1 + 2 * s / th2)^(-1.0) + (1 / th2) * (1 + 2 * s / th2)^(-1.5)) * exp(th2 * (1 - sqrt(1 + 2 * s / th2)))
    d2L <- w1 * d2L1 + w2 * d2L2
  } else if (frailty == "Gamma") {
    d2L <- ((theta + 1) / theta) * (1 + s / theta)^(-theta - 2)
  } else if (frailty == "IG") {
    d2L <- ((1 + 2 * s / theta)^(-1.0) + (1 / theta) * (1 + 2 * s / theta)^(-1.5)) * exp(theta * (1 - sqrt(1 + 2 * s / theta)))
  } else if (frailty == "Lindley") {
    d2L <- (2 * theta^2 * (s + theta + 3)) / ((1 + theta) * (s + theta)^4)
  } else if (frailty == "None") {
    d2L <- exp(-s)
  }
  return(pmax(d2L, 1e-15))
}

#' Variance of Frailty Distribution
#'
#' Computes the theoretical variance of the frailty distribution under the identifiability condition \eqn{E[W] = 1}.
#'
#' @param frailty Character string indicating frailty distribution.
#' @param theta Numeric parameter for frailty distribution.
#' @param mu Optional secondary parameter for \code{"GL"}.
#' @param theta2 Optional second parameter for \code{"GIG"}.
#'
#' @return Numeric scalar representing variance of frailty variable.
#' @export
#'
#' @examples
#' var_frailty("WL", theta = 3.4)
#' var_frailty("Gamma", theta = 2.0)
var_frailty <- function(frailty = c("WL", "GL", "GIG", "Gamma", "IG", "Lindley", "None"), theta, mu = NULL, theta2 = NULL) {
  frailty <- match.arg(frailty)
  
  if (frailty == "WL") {
    v <- 2 / (theta + sqrt(theta * (theta + 1)))
  } else if (frailty == "GL") {
    if (is.null(mu)) mu <- 0.5 * (1 + theta)
    num <- theta^4 - (theta^3) * mu + 3 * (theta^2) * (1 + theta) - 4 * (theta^2) * mu + 3 * theta * mu * (mu - 1) + mu^2
    den <- theta * (1 + theta)^2
    v <- num / den
  } else if (frailty == "GIG") {
    th2 <- if (is.null(theta2)) theta else theta2
    v <- (theta^3 + th2^3) / (theta + th2)
  } else if (frailty == "Gamma") {
    v <- 1 / theta
  } else if (frailty == "IG") {
    v <- 1 / theta
  } else if (frailty == "Lindley") {
    v <- (theta + 2) / (theta * (theta + 1))
  } else if (frailty == "None") {
    v <- 0.0
  }
  return(max(v, 0.0))
}

#' Kendall's Tau Measure of Dependence
#'
#' Computes Kendall's rank correlation coefficient \eqn{\tau} for the shared frailty model
#' under reversed hazard rate via numerical integration:
#' \deqn{\tau = 4 \int_0^\infty s L_W''(s) L_W(s) ds - 1}
#'
#' @param frailty Character string indicating frailty distribution.
#' @param theta Numeric parameter for frailty distribution.
#' @param mu Optional secondary parameter for \code{"GL"}.
#' @param theta2 Optional second parameter for \code{"GIG"}.
#'
#' @return Numeric scalar representing Kendall's tau in [0, 1].
#' @export
#'
#' @examples
#' kendall_tau_rhr("WL", theta = 3.4)
#' kendall_tau_rhr("Gamma", theta = 2.0)
kendall_tau_rhr <- function(frailty = c("WL", "GL", "GIG", "Gamma", "IG", "Lindley", "None"), theta, mu = NULL, theta2 = NULL) {
  frailty <- match.arg(frailty)
  if (frailty == "None") return(0.0)
  if (frailty == "Gamma") return(1 / (1 + 2 * theta))
  
  integrand <- function(s) {
    s * d2laplace_frailty(s, frailty, theta, mu, theta2) * laplace_frailty(s, frailty, theta, mu, theta2)
  }
  
  res <- try(stats::integrate(integrand, lower = 0, upper = Inf, rel.tol = 1e-6), silent = TRUE)
  if (inherits(res, "try-error")) {
    res <- stats::integrate(integrand, lower = 0, upper = 500, rel.tol = 1e-5)
  }
  tau <- 4 * res$value - 1
  return(min(max(tau, 0.0), 1.0))
}

#' Clayton-Oakes Cross-Ratio Function under Reversed Hazard Rate
#'
#' Computes the bivariate Clayton-Oakes cross-ratio association measure \eqn{\kappa(t_1, t_2)}.
#'
#' @param t1 Numeric scalar or vector of first component times.
#' @param t2 Numeric scalar or vector of second component times.
#' @param x1 Numeric matrix or vector of covariates for first component.
#' @param x2 Numeric matrix or vector of covariates for second component.
#' @param beta Numeric vector of regression coefficients.
#' @param frailty Character string indicating frailty distribution.
#' @param theta Numeric parameter for frailty distribution.
#' @param mu Optional secondary parameter for \code{"GL"}.
#' @param theta2 Optional second parameter for \code{"GIG"}.
#' @param baseline Character string indicating baseline distribution.
#' @param params1 Numeric named vector of baseline parameters for first component.
#' @param params2 Numeric named vector of baseline parameters for second component.
#'
#' @return Numeric vector of cross-ratio values.
#' @export
#'
#' @examples
#' cross_ratio_rhr(15, 18, frailty = "WL", theta = 3.4, baseline = "GE",
#'                 params1 = c(zeta = 12.4, delta = 0.11),
#'                 params2 = c(zeta = 11.8, delta = 0.11))
cross_ratio_rhr <- function(t1, t2, x1 = NULL, x2 = NULL, beta = NULL, frailty = "WL", theta = 1, mu = NULL, theta2 = NULL, baseline = "GE", params1 = c(zeta = 1, delta = 1), params2 = c(zeta = 1, delta = 1)) {
  R1 <- crhr_baseline(t1, baseline, params1)
  R2 <- crhr_baseline(t2, baseline, params2)
  
  eta1 <- if (!is.null(x1) && !is.null(beta)) exp(as.numeric(x1 %*% beta)) else rep(1, length(t1))
  eta2 <- if (!is.null(x2) && !is.null(beta)) exp(as.numeric(x2 %*% beta)) else rep(1, length(t2))
  
  s <- R1 * eta1 + R2 * eta2
  
  L <- laplace_frailty(s, frailty, theta, mu, theta2)
  dL <- dlaplace_frailty(s, frailty, theta, mu, theta2)
  d2L <- d2laplace_frailty(s, frailty, theta, mu, theta2)
  
  cr <- (d2L * L) / (dL^2)
  return(cr)
}

#' Conditional Probability Measure under Reversed Hazard Rate
#'
#' Computes the bivariate conditional dependence probability ratio \eqn{\chi(t_1, t_2) = S(t_1, t_2) / (S(t_1) S(t_2))}.
#'
#' @param t1 Numeric scalar or vector of first component times.
#' @param t2 Numeric scalar or vector of second component times.
#' @param x1 Numeric matrix or vector of covariates for first component.
#' @param x2 Numeric matrix or vector of covariates for second component.
#' @param beta Numeric vector of regression coefficients.
#' @param frailty Character string indicating frailty distribution.
#' @param theta Numeric parameter for frailty distribution.
#' @param mu Optional secondary parameter for \code{"GL"}.
#' @param theta2 Optional second parameter for \code{"GIG"}.
#' @param baseline Character string indicating baseline distribution.
#' @param params1 Numeric named vector of baseline parameters for first component.
#' @param params2 Numeric named vector of baseline parameters for second component.
#'
#' @return Numeric vector of conditional probability ratios.
#' @export
#'
#' @examples
#' cond_prob_rhr(15, 18, frailty = "WL", theta = 3.4, baseline = "GE",
#'               params1 = c(zeta = 12.4, delta = 0.11),
#'               params2 = c(zeta = 11.8, delta = 0.11))
cond_prob_rhr <- function(t1, t2, x1 = NULL, x2 = NULL, beta = NULL, frailty = "WL", theta = 1, mu = NULL, theta2 = NULL, baseline = "GE", params1 = c(zeta = 1, delta = 1), params2 = c(zeta = 1, delta = 1)) {
  R1 <- crhr_baseline(t1, baseline, params1)
  R2 <- crhr_baseline(t2, baseline, params2)
  
  eta1 <- if (!is.null(x1) && !is.null(beta)) exp(as.numeric(x1 %*% beta)) else rep(1, length(t1))
  eta2 <- if (!is.null(x2) && !is.null(beta)) exp(as.numeric(x2 %*% beta)) else rep(1, length(t2))
  
  s_joint <- R1 * eta1 + R2 * eta2
  s1 <- R1 * eta1
  s2 <- R2 * eta2
  
  F_joint <- laplace_frailty(s_joint, frailty, theta, mu, theta2)
  F1 <- laplace_frailty(s1, frailty, theta, mu, theta2)
  F2 <- laplace_frailty(s2, frailty, theta, mu, theta2)
  
  S_joint <- pmax(1 - F1 - F2 + F_joint, 1e-8)
  S1 <- pmax(1 - F1, 1e-8)
  S2 <- pmax(1 - F2, 1e-8)
  
  chi_val <- S_joint / (S1 * S2)
  return(chi_val)
}
