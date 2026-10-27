#' Dependence Measures and Theoretical Frailty Variance
#'
#' Computes Kendall's tau (rank-based association) and theoretical frailty
#' variance for additive shared frailty models under Generalized Lindley,
#' Weighted Lindley, Gamma, and Inverse Gaussian frailty distributions.
#'
#' @param frailty_type Character string specifying the frailty distribution:
#'   \code{"genlindley"}, \code{"wlindley"}, \code{"gamma"}, or \code{"invgaussian"}.
#' @param theta Positive parameter for Generalized Lindley, Gamma, or Inverse Gaussian frailty.
#' @param mu Positive parameter for Generalized Lindley frailty (with \code{0 < mu < 1 + theta}).
#' @param s_param Positive shape parameter for Weighted Lindley frailty.
#'
#' @return A scalar numeric value representing Kendall's tau or theoretical frailty variance.
#'
#' @name dependence
#' @references
#' Pandey, A., Hanagal, D. D., and Tyagi, S. (2022). Generalised Lindley shared
#' additive frailty regression model for bivariate survival data.
#' \emph{Statistics in Transition new series}, 23(4), 161-176.
#' \doi{10.2478/stattrans-2022-0048}
#'
#' Tyagi, S., Pandey, A., Hanagal, D. D., and Chesneau, C. (2024). A study on
#' comparisons of additive regression frailty models to counter heterogeneity:
#' Bayesian strategies and case study. \emph{Communications in Statistics -
#' Simulation and Computation}. \doi{10.1080/03610918.2024.2391874}
#'
#' @examples
#' # Kendall's tau for Generalized Lindley frailty
#' kendall_tau(frailty_type = "genlindley", theta = 2.29258, mu = 1.38391)
#'
#' # Kendall's tau for Weighted Lindley frailty
#' kendall_tau(frailty_type = "wlindley", s_param = 1.0955)
#'
#' # Frailty variances
#' frailty_variance(frailty_type = "genlindley", theta = 2.29, mu = 1.38)
#' frailty_variance(frailty_type = "wlindley", s_param = 1.10)
NULL

#' @rdname dependence
#' @export
kendall_tau <- function(frailty_type = c("genlindley", "wlindley", "gamma", "invgaussian"),
                        theta = NULL, mu = NULL, s_param = NULL) {
  frailty_type <- match.arg(frailty_type)

  if (frailty_type == "gamma") {
    if (is.null(theta) || theta <= 0) return(NA_real_)
    return(theta / (theta + 2))
  }

  if (frailty_type == "genlindley") {
    if (is.null(theta) || is.null(mu) || theta <= 0 || mu <= 0 || mu >= 1 + theta) {
      return(NA_real_)
    }
  } else if (frailty_type == "wlindley") {
    if (is.null(s_param) || s_param <= 0) {
      return(NA_real_)
    }
  } else if (frailty_type == "invgaussian") {
    if (is.null(theta) || theta <= 0) {
      return(NA_real_)
    }
  }

  integrand <- function(s) {
    s_pos <- pmax(s, 0)
    if (frailty_type == "genlindley") {
      L0 <- laplace_genlindley(s_pos, theta, mu)
      L2 <- laplace_deriv2_genlindley(s_pos, theta, mu)
    } else if (frailty_type == "wlindley") {
      L0 <- laplace_wlindley(s_pos, s_param)
      L2 <- laplace_deriv2_wlindley(s_pos, s_param)
    } else if (frailty_type == "invgaussian") {
      L0 <- laplace_invgaussian(s_pos, theta)
      L2 <- laplace_deriv2_invgaussian(s_pos, theta)
    }
    val <- 4 * s_pos * L2 * L0
    val[!is.finite(val)] <- 0
    val
  }

  res <- tryCatch({
    val <- stats::integrate(integrand, lower = 0, upper = Inf,
                            subdivisions = 500, rel.tol = 1e-6, stop.on.error = FALSE)
    if (!is.null(val$value) && is.finite(val$value)) {
      val$value - 1
    } else {
      # Fallback numerical integration on truncated interval [0, 500]
      val_trunc <- stats::integrate(integrand, lower = 0, upper = 500,
                                   subdivisions = 500, rel.tol = 1e-5, stop.on.error = FALSE)
      val_trunc$value - 1
    }
  }, error = function(e) {
    # Fallback to bounded domain integration
    val_b <- tryCatch({
      stats::integrate(integrand, lower = 0, upper = 200, rel.tol = 1e-4)$value - 1
    }, error = function(e2) NA_real_)
    val_b
  })

  if (is.na(res) || !is.finite(res)) return(NA_real_)
  pmin(pmax(res, 0), 1)
}

#' @rdname dependence
#' @export
frailty_variance <- function(frailty_type = c("genlindley", "wlindley", "gamma", "invgaussian"),
                             theta = NULL, mu = NULL, s_param = NULL) {
  frailty_type <- match.arg(frailty_type)
  if (frailty_type == "genlindley") {
    if (is.null(theta) || is.null(mu)) return(NA_real_)
    return(var_genlindley(theta, mu))
  } else if (frailty_type == "wlindley") {
    if (is.null(s_param)) return(NA_real_)
    return(var_wlindley(s_param))
  } else if (frailty_type == "gamma") {
    if (is.null(theta)) return(NA_real_)
    return(var_gamma(theta))
  } else if (frailty_type == "invgaussian") {
    if (is.null(theta)) return(NA_real_)
    return(var_invgaussian(theta))
  }
  NA_real_
}
