#' Bayesian Point Estimation and Posterior Risks under Multiple Loss Functions
#'
#' Computes Bayes point estimates and corresponding posterior risks under
#' Squared Error Loss Function (SELF), Modified Quadratic Squared Error
#' Loss Function (MQSELF), and Precautionary Loss Function (PLF), along with
#' equal-tailed Bayesian Credible Intervals and Highest Posterior Density (HPD) bounds.
#'
#' @param samples Matrix or data frame of MCMC posterior draws (rows = iterations, columns = parameters).
#' @param level Confidence level for credible intervals (default is 0.95).
#'
#' @return A data frame containing:
#' \describe{
#'   \item{Parameter}{Parameter name.}
#'   \item{Est_SELF}{Bayes estimate under SELF (posterior mean).}
#'   \item{Risk_SELF}{Posterior risk under SELF (posterior variance).}
#'   \item{Est_MQSELF}{Bayes estimate under MQSELF.}
#'   \item{Risk_MQSELF}{Posterior risk under MQSELF.}
#'   \item{Est_PLF}{Bayes estimate under PLF.}
#'   \item{Risk_PLF}{Posterior risk under PLF.}
#'   \item{LCL}{Lower bound of equal-tailed credible interval.}
#'   \item{UCL}{Upper bound of equal-tailed credible interval.}
#' }
#'
#' @export
#' @references
#' Tyagi, S., Pandey, A., Hanagal, D. D., and Chesneau, C. (2024). A study on
#' comparisons of additive regression frailty models to counter heterogeneity:
#' Bayesian strategies and case study. \emph{Communications in Statistics -
#' Simulation and Computation}. \doi{10.1080/03610918.2024.2391874}
#'
#' Norstrom, J. G. (1996). The use of precautionary loss functions in risk
#' analysis. \emph{IEEE Transactions on Reliability}, 45(3), 400-403.
#' \doi{10.1109/24.536992}
#'
#' @examples
#' # Simulated posterior samples
#' set.seed(123)
#' post_draws <- matrix(rgamma(1000 * 3, shape = 4, rate = 2), ncol = 3)
#' colnames(post_draws) <- c("zeta1", "delta1", "xi1")
#' calc_loss_estimates(post_draws)
calc_loss_estimates <- function(samples, level = 0.95) {
  samples <- as.matrix(samples)
  n_pars <- ncol(samples)
  par_names <- colnames(samples)
  if (is.null(par_names)) par_names <- paste0("param", seq_len(n_pars))

  alpha <- 1 - level
  prob_low <- alpha / 2
  prob_high <- 1 - alpha / 2

  res_list <- vector("list", n_pars)

  for (j in seq_len(n_pars)) {
    vals <- samples[, j]
    vals_clean <- vals[!is.na(vals) & !is.nan(vals)]
    if (length(vals_clean) == 0) vals_clean <- c(0, 0)

    # 1. SELF
    est_self <- mean(vals_clean)
    risk_self <- stats::var(vals_clean)

    # 2. MQSELF: E(1/theta) / E(1/theta^2)
    # Handle positive/negative parameters robustly
    if (all(vals_clean > 0)) {
      e_inv1 <- mean(1 / vals_clean)
      e_inv2 <- mean(1 / (vals_clean^2))
      est_mqself <- if (e_inv2 > 0) e_inv1 / e_inv2 else est_self
      risk_mqself <- if (e_inv2 > 0) pmax(0, 1 - (e_inv1^2) / e_inv2) else 0
    } else if (all(vals_clean < 0)) {
      abs_vals <- abs(vals_clean)
      e_inv1 <- mean(1 / abs_vals)
      e_inv2 <- mean(1 / (abs_vals^2))
      est_mqself <- if (e_inv2 > 0) - (e_inv1 / e_inv2) else est_self
      risk_mqself <- if (e_inv2 > 0) pmax(0, 1 - (e_inv1^2) / e_inv2) else 0
    } else {
      # Values span across zero: protect against division by zero
      shift <- max(abs(vals_clean)) + 1
      shifted <- vals_clean + shift
      e_inv1 <- mean(1 / shifted)
      e_inv2 <- mean(1 / (shifted^2))
      est_mqself <- (e_inv1 / e_inv2) - shift
      risk_mqself <- pmax(0, 1 - (e_inv1^2) / e_inv2)
    }

    # 3. PLF: sqrt(E(theta^2))
    e_sq <- mean(vals_clean^2)
    est_plf <- sqrt(pmax(0, e_sq))
    # If parameter is predominantly negative, preserve direction
    if (est_self < 0) {
      est_plf_signed <- -est_plf
      risk_plf <- 2 * abs(est_plf - abs(est_self))
    } else {
      est_plf_signed <- est_plf
      risk_plf <- 2 * pmax(0, est_plf - est_self)
    }

    # Credible Intervals
    ci_bounds <- stats::quantile(vals_clean, probs = c(prob_low, prob_high), names = FALSE)

    res_list[[j]] <- data.frame(
      Parameter = par_names[j],
      Est_SELF = est_self,
      Risk_SELF = risk_self,
      Est_MQSELF = est_mqself,
      Risk_MQSELF = risk_mqself,
      Est_PLF = est_plf_signed,
      Risk_PLF = risk_plf,
      LCL = ci_bounds[1],
      UCL = ci_bounds[2],
      stringsAsFactors = FALSE
    )
  }

  do.call(rbind, res_list)
}
