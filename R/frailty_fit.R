#' Unified Formula Interface for Frailty Models
#'
#' Fits frailty survival regression models using either Bayesian MCMC or Maximum
#' Likelihood Estimation.
#'
#' @param formula Formula expression for bivariate or univariate survival data.
#' @param data Data frame.
#' @param engine Estimation engine: \code{"bayes"} (default) or \code{"mle"}.
#' @param ... Arguments passed to \code{\link{frailty_bayes}} or \code{\link{frailty_mle}}.
#'
#' @return A fitted model object of class \code{"frailty_bayes"} or \code{"frailty_mle"}.
#' @export
#' @examples
#' data(kidney_frailty)
#' fit <- frailty_fit(
#'   Surv(time1, time2, status1, status2) ~ age + sex,
#'   data = kidney_frailty,
#'   engine = "mle",
#'   baseline = "weibull",
#'   frailty = "gamma"
#' )
#' print(fit)
frailty_fit <- function(formula, data, engine = c("bayes", "mle"), ...) {
  engine <- match.arg(engine)
  if (engine == "bayes") {
    frailty_bayes(formula, data, ...)
  } else {
    frailty_mle(formula, data, ...)
  }
}
