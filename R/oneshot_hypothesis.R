#' Likelihood Ratio Test for Component Independence
#'
#' @description
#' Performs the Likelihood Ratio Test (LRT) to evaluate the null hypothesis
#' of component independence against the shared gamma frailty alternative.
#'
#' @param object A fitted model object of class \code{"oneshot_kofm"} or \code{"oneshot_frailty"}.
#' @param ... Additional arguments passed to methods (such as \code{fit_indep} for \code{"oneshot_frailty"}).
#'
#' @return An S3 object of class \code{"htest"} containing the test statistic, degrees of freedom, and p-value.
#' @export
#' @examples
#' data(classH)
#' fit_w <- oneshot_kofm(
#'   data = classH,
#'   time = "time",
#'   stress = "stress",
#'   counts = c("n_none", "n_T", "n_G", "n_both"),
#'   baseline = "exponential",
#'   frailty = TRUE
#' )
#' test_independence(fit_w)
test_independence <- function(object, ...) {
  UseMethod("test_independence")
}

#' @export
test_independence.oneshot_kofm <- function(object, ...) {
  kofm_test_independence(object = object)
}

#' @export
test_independence.oneshot_frailty <- function(object, fit_indep = NULL, ...) {
  frailty_test_independence(fit_frailty = object, fit_indep = fit_indep)
}

#' Likelihood Ratio Test for Baseline Lifetime Distribution
#'
#' @description
#' Performs the Likelihood Ratio Test comparing Weibull versus Exponential baseline lifetime distributions.
#'
#' @param object A fitted Weibull model object of class \code{"oneshot_kofm"} or \code{"oneshot_frailty"}.
#' @param ... Additional arguments passed to methods (such as \code{fit_exp} for \code{"oneshot_frailty"}).
#'
#' @return An S3 object of class \code{"htest"} containing the test statistic, degrees of freedom, and p-value.
#' @export
#' @examples
#' data(classH)
#' fit_w <- oneshot_kofm(
#'   data = classH,
#'   time = "time",
#'   stress = "stress",
#'   counts = c("n_none", "n_T", "n_G", "n_both"),
#'   baseline = "weibull",
#'   frailty = TRUE
#' )
#' test_baseline(fit_w)
test_baseline <- function(object, ...) {
  UseMethod("test_baseline")
}

#' @export
test_baseline.oneshot_kofm <- function(object, ...) {
  kofm_test_baseline(object = object)
}

#' @export
test_baseline.oneshot_frailty <- function(object, fit_exp = NULL, ...) {
  frailty_test_baseline(fit_weibull = object, fit_exponential = fit_exp)
}


#' Likelihood Ratio Test for Component Independence
#'
#' @description
#' Performs the Likelihood Ratio Test (LRT) and Wald test for testing the null hypothesis
#' of component independence (\eqn{H_0: \beta = 0}) against the alternative of positive
#' dependence under shared gamma frailty (\eqn{H_1: \beta > 0}) as derived in Ling et al. (2021).
#'
#' @param object A fitted model object of class \code{"oneshot_kofm"} with \code{frailty = TRUE}.
#'
#' @return An S3 object of class \code{"htest"} containing:
#'   \describe{
#'     \item{statistic}{The LRT test statistic \eqn{\Lambda_{\text{LRT}} = -2(\ell_{\text{ind}} - \ell_{\text{frailty}})}.}
#'     \item{p.value}{The asymptotic \eqn{p}-value based on the 50:50 mixture chi-square null distribution \eqn{0.5 \chi_0^2 + 0.5 \chi_1^2}.}
#'     \item{parameter}{Degrees of freedom (\code{df = 1}).}
#'     \item{method}{Character string describing the test.}
#'     \item{data.name}{Character string naming the input dataset.}
#'     \item{estimate}{Estimated frailty variance \eqn{\hat{\beta}}.}
#'     \item{null.value}{The hypothesized value \eqn{\beta = 0}.}
#'     \item{alternative}{Direction of alternative hypothesis (\code{"greater"}).}
#'     \item{wald}{A list containing the Wald test statistic, standard error, and one-sided \eqn{p}-value.}
#'   }
#'
#' @details
#' Since \eqn{\beta = 0} lies on the boundary of the parameter space, the asymptotic null
#' distribution of the LRT statistic is a mixture of chi-square distributions:
#' \deqn{P(\Lambda_{\text{LRT}} \le \lambda) = 0.5 + 0.5 P(\chi_1^2 \le \lambda)}
#' and the \eqn{p}-value is:
#' \deqn{p = 0.5 \times P(\chi_1^2 \ge \Lambda_{\text{LRT}})}
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
#' fit <- oneshot_kofm(
#'   data = classH,
#'   time = "time",
#'   stress = "stress",
#'   counts = c("n_none", "n_T", "n_G", "n_both"),
#'   baseline = "exponential",
#'   frailty = TRUE
#' )
#' test_independence(fit)
#' @noRd
kofm_test_independence <- function(object) {
  if (!inherits(object, "oneshot_kofm")) {
    stop("'object' must be a fitted object of class 'oneshot_kofm'.")
  }
  
  if (!object$frailty) {
    stop("The supplied model was fitted with 'frailty = FALSE'. Fit a model with 'frailty = TRUE' to test independence.")
  }
  
  # Fit independence model
  call_ind <- object$call
  call_ind$frailty <- FALSE
  call_ind$start <- NULL
  
  fit_ind <- tryCatch({
    eval(call_ind, parent.frame())
  }, error = function(e) {
    oneshot_kofm(
      data = object$data,
      time = as.character(call_ind$time),
      stress = if (!is.null(call_ind$stress)) as.character(call_ind$stress) else NULL,
      counts = if (!is.null(call_ind$counts)) as.character(call_ind$counts) else NULL,
      components = object$components,
      baseline = object$baseline,
      frailty = FALSE,
      method = object$method
    )
  })
  
  ll_frailty <- object$loglik
  ll_ind <- fit_ind$loglik
  
  lrt_stat <- max(0, -2 * (ll_ind - ll_frailty))
  
  # 50:50 mixture chi-square p-value
  p_val_lrt <- if (lrt_stat > 0) 0.5 * stats::pchisq(lrt_stat, df = 1, lower.tail = FALSE) else 1.0
  
  # Wald test
  beta_est <- object$coefficients["beta"]
  beta_se <- object$se["beta"]
  wald_stat <- if (beta_se > 0) as.numeric(beta_est / beta_se) else 0
  p_val_wald <- stats::pnorm(wald_stat, lower.tail = FALSE)
  
  data_name <- deparse(object$call$data)
  if (length(data_name) == 0) data_name <- "one-shot device data"
  
  structure(
    list(
      statistic = c("LRT" = lrt_stat),
      p.value = p_val_lrt,
      parameter = c("df" = 1),
      method = "Likelihood Ratio Test for Component Independence (Mixture Chi-Square)",
      data.name = data_name,
      estimate = c("beta" = beta_est),
      null.value = c("beta" = 0),
      alternative = "greater",
      loglik_frailty = ll_frailty,
      loglik_independence = ll_ind,
      wald = list(
        statistic = c("Wald z" = wald_stat),
        se = beta_se,
        p.value = p_val_wald
      )
    ),
    class = "htest"
  )
}

#' Likelihood Ratio Test for Baseline Distribution Discrimination
#'
#' @description
#' Performs a Likelihood Ratio Test comparing the Weibull baseline lifetime model
#' against the Exponential baseline model (\eqn{H_0: \alpha_m = 1} for all \eqn{m}).
#'
#' @param object A fitted model object of class \code{"oneshot_kofm"} with \code{baseline = "weibull"}.
#'
#' @return An S3 object of class \code{"htest"} containing:
#'   \describe{
#'     \item{statistic}{The LRT test statistic \eqn{-2(\ell_{\text{exp}} - \ell_{\text{weibull}})}.}
#'     \item{p.value}{The asymptotic \eqn{p}-value from the \eqn{\chi^2} distribution.}
#'     \item{parameter}{Degrees of freedom (\eqn{M \times p}).}
#'     \item{method}{Character string describing the test.}
#'     \item{data.name}{Character string naming the input dataset.}
#'   }
#'
#' @details
#' Missing values are not permitted.
#'
#' @references
#' Ling, M.-H., Balakrishnan, N., Yu, C., and So, H. Y. (2021). Inference for One-Shot
#' Devices with Dependent k-Out-of-M Structured Components under Gamma Frailty.
#' \emph{Mathematics}, 9(23), 3032. \doi{10.3390/math9233032}.
#'
#' @examples
#' data(classH)
#' fit_w <- oneshot_kofm(
#'   data = classH,
#'   time = "time",
#'   stress = "stress",
#'   counts = c("n_none", "n_T", "n_G", "n_both"),
#'   baseline = "weibull",
#'   frailty = TRUE
#' )
#' test_baseline(fit_w)
#' @noRd
kofm_test_baseline <- function(object) {
  if (!inherits(object, "oneshot_kofm")) {
    stop("'object' must be a fitted object of class 'oneshot_kofm'.")
  }
  
  if (object$baseline != "weibull") {
    stop("The supplied model was fitted with 'baseline = \"exponential\"'. Supply a Weibull model to test baseline.")
  }
  
  call_exp <- object$call
  call_exp$baseline <- "exponential"
  call_exp$start <- NULL
  
  fit_exp <- tryCatch({
    eval(call_exp, parent.frame())
  }, error = function(e) {
    oneshot_kofm(
      data = object$data,
      time = as.character(call_exp$time),
      stress = if (!is.null(call_exp$stress)) as.character(call_exp$stress) else NULL,
      counts = if (!is.null(call_exp$counts)) as.character(call_exp$counts) else NULL,
      components = object$components,
      baseline = "exponential",
      frailty = object$frailty,
      method = object$method
    )
  })
  
  ll_w <- object$loglik
  ll_e <- fit_exp$loglik
  
  df_diff <- object$n_params - fit_exp$n_params
  lrt_stat <- max(0, -2 * (ll_e - ll_w))
  p_val <- stats::pchisq(lrt_stat, df = df_diff, lower.tail = FALSE)
  
  data_name <- deparse(object$call$data)
  if (length(data_name) == 0) data_name <- "one-shot device data"
  
  structure(
    list(
      statistic = c("LRT" = lrt_stat),
      p.value = p_val,
      parameter = c("df" = df_diff),
      method = "Likelihood Ratio Test for Baseline Distribution (Weibull vs. Exponential)",
      data.name = data_name,
      null.value = c("shape_parameters" = 1),
      alternative = "Weibull baseline fits significantly better"
    ),
    class = "htest"
  )
}


#' Likelihood Ratio Test for Component Independence
#'
#' @description
#' Performs the Likelihood Ratio Test (LRT) to evaluate the hypothesis of component
#' independence (\eqn{H_0: \beta = 0}) against the shared gamma frailty model (\eqn{H_1: \beta > 0}).
#'
#' @param fit_frailty An object of class \code{"oneshot_frailty"} fitted with \code{frailty = TRUE}.
#' @param fit_indep Optional object of class \code{"oneshot_frailty"} fitted with \code{frailty = FALSE}.
#'   If not provided, the independence model is fitted automatically.
#'
#' @return A list of class \code{"htest"} containing:
#'   \describe{
#'     \item{statistic}{The LRT test statistic \eqn{\Lambda = -2(\log L_{0} - \log L_{1})}.}
#'     \item{p.value}{The asymptotic p-value calculated from the mixture chi-square distribution \eqn{0.5 \chi_0^2 + 0.5 \chi_1^2}.}
#'     \item{estimate}{Estimated frailty variance \eqn{\beta} and Kendall's tau.}
#'     \item{method}{Character string describing the test.}
#'     \item{data.name}{Character string specifying the dataset name.}
#'     \item{null.value}{The null hypothesis value (\eqn{\beta = 0}).}
#'     \item{alternative}{Character string: \code{"greater"} (\eqn{\beta > 0}).}
#'   }
#'
#' @details
#' Because the null hypothesis \eqn{\beta = 0} lies on the boundary of the parameter space,
#' the standard regularity conditions do not apply. Yu, So, and Balakrishnan (2025) show
#' that the asymptotic distribution of the test statistic \eqn{\Lambda} under \eqn{H_0}
#' is a 50:50 mixture of a point mass at zero (\eqn{\chi_0^2}) and a 1-df chi-square (\eqn{\chi_1^2}):
#' \deqn{P(\Lambda \le \lambda) = 0.5 + 0.5 P(\chi_1^2 \le \lambda)}
#'
#' Missing values in input models are not permitted.
#'
#' @references
#' Yu, C., So, H. Y., and Balakrishnan, N. (2025). \emph{Inference for one-shot devices
#' with Weibull component lifetimes under gamma frailty model}. In: Stochastic Modeling
#' and Statistical Methods, Elsevier Inc., Chapter 1, pp. 1--34.
#' \doi{10.1016/B978-0-44-331694-4.00006-2}.
#'
#' @examples
#' data(classH)
#' fit_w <- oneshot_frailty(data = classH, time = "time", stress = "stress",
#'                          counts = c("n_none", "n_T", "n_G", "n_both"),
#'                          baseline = "weibull", frailty = TRUE)
#' test_independence(fit_w)
#' @noRd
frailty_test_independence <- function(fit_frailty, fit_indep = NULL) {
  if (!inherits(fit_frailty, "oneshot_frailty")) {
    stop("'fit_frailty' must be an object of class 'oneshot_frailty'.")
  }
  if (!fit_frailty$frailty) {
    stop("'fit_frailty' must have been fitted with 'frailty = TRUE'.")
  }
  
  if (is.null(fit_indep)) {
    fit_indep <- oneshot_frailty(
      data = fit_frailty$data,
      time = as.character(fit_frailty$call$time %||% "time"),
      stress = if (ncol(fit_frailty$X_cov) > 1) colnames(fit_frailty$X_cov)[-1] else NULL,
      counts = colnames(fit_frailty$count_mat),
      components = fit_frailty$components,
      baseline = fit_frailty$baseline,
      frailty = FALSE,
      method = fit_frailty$method
    )
  }
  
  ll_frailty <- fit_frailty$loglik
  ll_indep <- fit_indep$loglik
  
  lrt_stat <- max(0, -2 * (ll_indep - ll_frailty))
  names(lrt_stat) <- "LRT"
  
  if (lrt_stat > 0) {
    p_val <- 0.5 * stats::pchisq(lrt_stat, df = 1, lower.tail = FALSE)
  } else {
    p_val <- 1.0
  }
  
  beta_est <- fit_frailty$coefficients["beta"]
  tau_est <- fit_frailty$kendall_tau$tau
  est_vec <- c(beta = beta_est, kendall_tau = tau_est)
  
  structure(
    list(
      statistic = lrt_stat,
      p.value = p_val,
      estimate = est_vec,
      null.value = c(beta = 0),
      alternative = "greater (components are correlated via gamma frailty)",
      method = "Likelihood Ratio Test for Component Independence (Mixture Chi-Square)",
      data.name = deparse(substitute(fit_frailty))
    ),
    class = "htest"
  )
}

#' Likelihood Ratio Test for Baseline Hazard Distribution
#'
#' @description
#' Evaluates whether the component lifetime baseline hazard is Exponential (\eqn{H_0: \alpha_{i,m} = 1})
#' versus Weibull (\eqn{H_1: \alpha_{i,m} \ne 1}) using a Likelihood Ratio Test.
#'
#' @param fit_weibull An object of class \code{"oneshot_frailty"} fitted with \code{baseline = "weibull"}.
#' @param fit_exponential Optional object of class \code{"oneshot_frailty"} fitted with \code{baseline = "exponential"}.
#'
#' @return An object of class \code{"htest"} containing:
#'   \describe{
#'     \item{statistic}{The LRT test statistic \eqn{\Lambda = -2(\log L_{0} - \log L_{1})}.}
#'     \item{parameter}{Degrees of freedom for the chi-square test.}
#'     \item{p.value}{The p-value from the asymptotic \eqn{\chi^2} distribution.}
#'     \item{method}{Character string describing the test.}
#'     \item{data.name}{Name of the model object.}
#'   }
#'
#' @details
#' Under \eqn{H_0}, the shape regression parameters \eqn{b_{m,0} = b_{m,1} = \dots = 0}.
#' The asymptotic null distribution of \eqn{\Lambda} is \eqn{\chi^2_{df}} where \eqn{df}
#' equals the number of restricted shape parameters (\eqn{M x p_{cov}}).
#'
#' @references
#' Yu, C., So, H. Y., and Balakrishnan, N. (2025). \emph{Inference for one-shot devices
#' with Weibull component lifetimes under gamma frailty model}. In: Stochastic Modeling
#' and Statistical Methods, Elsevier Inc., Chapter 1, pp. 1--34.
#' \doi{10.1016/B978-0-44-331694-4.00006-2}.
#'
#' @examples
#' data(classH)
#' fit_w <- oneshot_frailty(data = classH, time = "time", stress = "stress",
#'                          counts = c("n_none", "n_T", "n_G", "n_both"),
#'                          baseline = "weibull", frailty = TRUE)
#' test_baseline(fit_w)
#' @noRd
frailty_test_baseline <- function(fit_weibull, fit_exponential = NULL) {
  if (!inherits(fit_weibull, "oneshot_frailty")) {
    stop("'fit_weibull' must be an object of class 'oneshot_frailty'.")
  }
  if (fit_weibull$baseline != "weibull") {
    stop("'fit_weibull' must have been fitted with baseline = 'weibull'.")
  }
  
  if (is.null(fit_exponential)) {
    fit_exponential <- oneshot_frailty(
      data = fit_weibull$data,
      time = as.character(fit_weibull$call$time %||% "time"),
      stress = if (ncol(fit_weibull$X_cov) > 1) colnames(fit_weibull$X_cov)[-1] else NULL,
      counts = colnames(fit_weibull$count_mat),
      components = fit_weibull$components,
      baseline = "exponential",
      frailty = fit_weibull$frailty,
      method = fit_weibull$method
    )
  }
  
  ll_weib <- fit_weibull$loglik
  ll_exp <- fit_exponential$loglik
  
  lrt_stat <- max(0, -2 * (ll_exp - ll_weib))
  df_val <- fit_weibull$n_params - fit_exponential$n_params
  names(lrt_stat) <- "LRT"
  names(df_val) <- "df"
  
  p_val <- stats::pchisq(lrt_stat, df = df_val, lower.tail = FALSE)
  
  structure(
    list(
      statistic = lrt_stat,
      parameter = df_val,
      p.value = p_val,
      method = "Likelihood Ratio Test for Baseline Hazard (Exponential vs Weibull)",
      data.name = deparse(substitute(fit_weibull))
    ),
    class = "htest"
  )
}
