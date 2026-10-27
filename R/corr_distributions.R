#' Baseline Hazard and Cumulative Hazard Functions
#'
#' @description
#' Internal and utility functions to evaluate baseline cumulative hazard H_0(t),
#' baseline hazard rate h_0(t), parameter names, and inverse cumulative hazard functions
#' for the 7 supported baseline lifetime distributions.
#'
#' @name distributions
#' @keywords internal
NULL

baseline_param_names <- function(dist) {
  dist <- match.arg(dist, c("genweibull", "genloglogistic", "weibull",
                            "exponential", "pareto", "lfr", "burr"))
  switch(dist,
    genweibull = c("alpha", "lambda", "gamma"),
    genloglogistic = c("alpha", "lambda", "gamma"),
    weibull = c("lambda", "gamma"),
    exponential = c("lambda"),
    pareto = c("alpha", "lambda"),
    lfr = c("alpha", "lambda"),
    burr = c("alpha", "lambda")
  )
}

baseline_param_count <- function(dist) {
  length(baseline_param_names(dist))
}

baseline_H_h <- function(t, dist, params) {
  t <- pmax(as.numeric(t), 1e-12)
  dist <- match.arg(dist, c("genweibull", "genloglogistic", "weibull",
                            "exponential", "pareto", "lfr", "burr"))
  
  if (dist == "genweibull") {
    alpha <- params[1]; lam <- params[2]; gamma <- params[3]
    alpha <- pmax(alpha, 1e-10); lam <- pmax(lam, 1e-10); gamma <- pmax(gamma, 1e-10)
    z <- exp(-lam * (t^gamma))
    u <- pmax(1e-12, 1 - z)
    u_alpha <- pmin(1 - 1e-12, u^alpha)
    H <- -log(pmax(1e-12, 1 - u_alpha))
    h <- (alpha * lam * gamma * (t^(gamma - 1)) * z * (u^(alpha - 1))) / pmax(1e-12, 1 - u_alpha)
  } else if (dist == "genloglogistic") {
    alpha <- params[1]; lam <- params[2]; gamma <- params[3]
    alpha <- pmax(alpha, 1e-10); lam <- pmax(lam, 1e-10); gamma <- pmax(gamma, 1e-10)
    arg <- 1 + lam * (t^gamma)
    H <- alpha * log(arg)
    h <- (alpha * lam * gamma * (t^(gamma - 1))) / arg
  } else if (dist == "weibull") {
    lam <- params[1]; gamma <- params[2]
    lam <- pmax(lam, 1e-10); gamma <- pmax(gamma, 1e-10)
    H <- lam * (t^gamma)
    h <- lam * gamma * (t^(gamma - 1))
  } else if (dist == "exponential") {
    lam <- params[1]
    lam <- pmax(lam, 1e-10)
    H <- lam * t
    h <- rep(lam, length(t))
  } else if (dist == "pareto") {
    alpha <- params[1]; lam <- params[2]
    alpha <- pmax(alpha, 1e-10); lam <- pmax(lam, 1e-10)
    H <- alpha * log(1 + lam * t)
    h <- (alpha * lam) / (1 + lam * t)
  } else if (dist == "lfr") {
    alpha <- params[1]; lam <- params[2]
    alpha <- pmax(alpha, 1e-10); lam <- pmax(lam, 1e-10)
    H <- alpha * t + 0.5 * lam * (t^2)
    h <- alpha + lam * t
  } else if (dist == "burr") {
    alpha <- params[1]; lam <- params[2]
    alpha <- pmax(alpha, 1e-10); lam <- pmax(lam, 1e-10)
    H <- alpha * log(1 + t^lam)
    h <- (alpha * lam * (t^(lam - 1))) / (1 + t^lam)
  }
  
  H <- pmax(H, 1e-12)
  h <- pmax(h, 1e-12)
  list(H = H, h = h)
}

baseline_inv_H <- function(u, dist, params) {
  u <- pmax(as.numeric(u), 1e-12)
  dist <- match.arg(dist, c("genweibull", "genloglogistic", "weibull",
                            "exponential", "pareto", "lfr", "burr"))
  
  if (dist == "genweibull") {
    alpha <- params[1]; lam <- params[2]; gamma <- params[3]
    val1 <- pmin(1 - 1e-12, 1 - exp(-u))
    val2 <- pmin(1 - 1e-12, val1^(1 / alpha))
    inner <- -log(pmax(1e-12, 1 - val2)) / lam
    t <- pmax(1e-12, inner)^(1 / gamma)
  } else if (dist == "genloglogistic") {
    alpha <- params[1]; lam <- params[2]; gamma <- params[3]
    inner <- pmax(0, exp(u / alpha) - 1) / lam
    t <- pmax(1e-12, inner)^(1 / gamma)
  } else if (dist == "weibull") {
    lam <- params[1]; gamma <- params[2]
    t <- (u / lam)^(1 / gamma)
  } else if (dist == "exponential") {
    lam <- params[1]
    t <- u / lam
  } else if (dist == "pareto") {
    alpha <- params[1]; lam <- params[2]
    t <- pmax(0, exp(u / alpha) - 1) / lam
  } else if (dist == "lfr") {
    alpha <- params[1]; lam <- params[2]
    t <- (-alpha + sqrt(pmax(0, alpha^2 + 2 * lam * u))) / lam
  } else if (dist == "burr") {
    alpha <- params[1]; lam <- params[2]
    inner <- pmax(0, exp(u / alpha) - 1)
    t <- pmax(1e-12, inner)^(1 / lam)
  }
  pmax(t, 1e-12)
}
