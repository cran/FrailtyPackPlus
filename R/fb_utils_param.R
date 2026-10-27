#' Parse Data and Formulas for Frailty Models
#' @param formula Formula expression.
#' @param data Data frame.
#' @param type Model type string.
#' @return Parsed data list.
#' @keywords internal
parse_frailty_data <- function(formula, data, type) {
  if (is.null(data) || !is.data.frame(data)) {
    stop("'data' must be a non-null data frame.")
  }
  
  f_terms <- stats::terms(formula, data = data)
  lhs <- formula[[2]]
  rhs_formula <- stats::as.formula(paste("~", paste(attr(f_terms, "term.labels"), collapse = " + ")))
  
  # Covariate model matrix
  if (length(attr(f_terms, "term.labels")) > 0) {
    cov_mat <- stats::model.matrix(rhs_formula, data = data)
    if (ncol(cov_mat) > 1 && colnames(cov_mat)[1] == "(Intercept)") {
      cov_mat <- cov_mat[, -1, drop = FALSE]
    } else if (ncol(cov_mat) == 1 && colnames(cov_mat)[1] == "(Intercept)") {
      cov_mat <- matrix(0, nrow = nrow(data), ncol = 0)
    }
  } else {
    cov_mat <- matrix(0, nrow = nrow(data), ncol = 0)
  }
  
  p_cov <- ncol(cov_mat)
  covar_names <- colnames(cov_mat)
  n <- nrow(data)
  
  # Parse LHS elements
  lhs_list <- as.list(lhs)
  fn_name <- as.character(lhs_list[[1]])
  
  if (grepl("bivariate", type)) {
    if (length(lhs_list) >= 5) {
      # Surv(time1, time2, status1, status2) or cbind(time1, time2, status1, status2)
      t1 <- eval(lhs_list[[2]], envir = data, enclos = parent.frame())
      t2 <- eval(lhs_list[[3]], envir = data, enclos = parent.frame())
      d1 <- eval(lhs_list[[4]], envir = data, enclos = parent.frame())
      d2 <- eval(lhs_list[[5]], envir = data, enclos = parent.frame())
    } else if (all(c("time1", "time2", "status1", "status2") %in% names(data))) {
      t1 <- data[["time1"]]
      t2 <- data[["time2"]]
      d1 <- data[["status1"]]
      d2 <- data[["status2"]]
    } else {
      stop("For bivariate models, provide Surv(time1, time2, status1, status2) or include columns time1, time2, status1, status2 in data.")
    }
    
    list(n = n, time1 = as.numeric(t1), time2 = as.numeric(t2),
         status1 = as.numeric(d1), status2 = as.numeric(d2),
         X1 = cov_mat, X2 = cov_mat, p_cov = p_cov, covar_names = covar_names)
  } else {
    # Univariate
    if (length(lhs_list) >= 3) {
      t <- eval(lhs_list[[2]], envir = data, enclos = parent.frame())
      d <- eval(lhs_list[[3]], envir = data, enclos = parent.frame())
    } else if (all(c("time", "status") %in% names(data))) {
      t <- data[["time"]]
      d <- data[["status"]]
    } else {
      stop("For univariate models, provide Surv(time, status) or include columns time, status in data.")
    }
    
    list(n = n, time = as.numeric(t), status = as.numeric(d),
         X = cov_mat, p_cov = p_cov, covar_names = covar_names)
  }
}

#' Parameter Names Generator
#' @param baseline Baseline name.
#' @param frailty Frailty name.
#' @param covar_names Covariate names.
#' @return Character vector of parameter names.
#' @keywords internal
get_param_names <- function(baseline, frailty, covar_names) {
  base_names <- switch(baseline,
    weibull = c("lambda", "rho"),
    exponential = c("lambda"),
    genexponential = c("lambda", "alpha"),
    exppower = c("lambda", "alpha"),
    lfr = c("a", "b"),
    genloglogistic = c("lambda", "alpha", "gamma"),
    genweibull = c("lambda", "alpha", "gamma"),
    genrayleigh = c("lambda", "alpha"),
    weightedexp = c("lambda", "alpha"),
    extendedweibull = c("lambda", "rho")
  )
  
  frail_names <- switch(frailty,
    gamma = c("theta"),
    invgauss = c("theta"),
    posstable = c("alpha_frail"),
    compoundpoisson = c("theta1", "theta2"),
    compoundnb = c("theta1", "theta2", "alpha_nb"),
    genlindley = c("theta", "alpha_lindley"),
    weightedlindley = c("theta", "c_param"),
    mixture_invgauss = c("p1", "theta1", "theta2")
  )
  
  c(paste0("base_", base_names), covar_names, paste0("frail_", frail_names))
}

#' Parameter Transformations for Constrained to Unconstrained Mapping
#' @param baseline Baseline name.
#' @param frailty Frailty name.
#' @param n_cov Number of covariates.
#' @return List of transformation functions.
#' @keywords internal
get_param_transforms <- function(baseline, frailty, n_cov) {
  n_base <- get_n_base_params(baseline)
  n_frail <- get_n_frail_params(frailty)
  
  to_uncon <- function(con_p) {
    uncon <- numeric(length(con_p))
    uncon[1:n_base] <- log(pmax(con_p[1:n_base], 1e-8))
    if (n_cov > 0) {
      uncon[(n_base + 1):(n_base + n_cov)] <- con_p[(n_base + 1):(n_base + n_cov)]
    }
    idx_f <- (n_base + n_cov + 1):(n_base + n_cov + n_frail)
    if (frailty == "posstable") {
      a <- min(max(con_p[idx_f[1]], 1e-6), 0.9999)
      uncon[idx_f[1]] <- stats::qlogis(a)
    } else if (frailty == "mixture_invgauss") {
      p1 <- min(max(con_p[idx_f[1]], 1e-6), 0.9999)
      uncon[idx_f[1]] <- stats::qlogis(p1)
      uncon[idx_f[2]] <- log(max(con_p[idx_f[2]], 1e-8))
      uncon[idx_f[3]] <- log(max(con_p[idx_f[3]], 1e-8))
    } else {
      uncon[idx_f] <- log(pmax(con_p[idx_f], 1e-8))
    }
    uncon
  }
  
  to_con <- function(uncon_p) {
    con <- numeric(length(uncon_p))
    con[1:n_base] <- exp(pmin(pmax(uncon_p[1:n_base], -20), 20))
    if (n_cov > 0) {
      con[(n_base + 1):(n_base + n_cov)] <- uncon_p[(n_base + 1):(n_base + n_cov)]
    }
    idx_f <- (n_base + n_cov + 1):(n_base + n_cov + n_frail)
    if (frailty == "posstable") {
      con[idx_f[1]] <- stats::plogis(uncon_p[idx_f[1]])
    } else if (frailty == "mixture_invgauss") {
      con[idx_f[1]] <- stats::plogis(uncon_p[idx_f[1]])
      con[idx_f[2]] <- exp(pmin(pmax(uncon_p[idx_f[2]], -20), 20))
      con[idx_f[3]] <- exp(pmin(pmax(uncon_p[idx_f[3]], -20), 20))
    } else {
      con[idx_f] <- exp(pmin(pmax(uncon_p[idx_f], -20), 20))
    }
    con
  }
  
  log_jacobian <- function(con_p) {
    jac <- sum(log(pmax(con_p[1:n_base], 1e-8)))
    idx_f <- (n_base + n_cov + 1):(n_base + n_cov + n_frail)
    if (frailty == "posstable") {
      a <- con_p[idx_f[1]]
      jac <- jac + log(pmax(a * (1 - a), 1e-12))
    } else if (frailty == "mixture_invgauss") {
      p1 <- con_p[idx_f[1]]
      jac <- jac + log(pmax(p1 * (1 - p1), 1e-12)) + log(pmax(con_p[idx_f[2]], 1e-8)) + log(pmax(con_p[idx_f[3]], 1e-8))
    } else {
      jac <- jac + sum(log(pmax(con_p[idx_f], 1e-8)))
    }
    jac
  }
  
  list(to_uncon = to_uncon, to_con = to_con, log_jacobian = log_jacobian)
}

#' Prior Distribution Evaluator
#' @param con_p Constrained parameters.
#' @param n_base Number of baseline parameters.
#' @param n_cov Number of covariates.
#' @param n_frail Number of frailty parameters.
#' @param baseline Baseline name.
#' @param frailty Frailty name.
#' @param priors Optional priors list.
#' @return Numeric scalar log-prior.
#' @keywords internal
eval_log_prior <- function(con_p, n_base, n_cov, n_frail, baseline, frailty, priors = NULL) {
  lp <- 0
  for (i in 1:n_base) {
    lp <- lp + stats::dgamma(con_p[i], shape = 0.01, rate = 0.01, log = TRUE)
  }
  if (n_cov > 0) {
    for (j in 1:n_cov) {
      idx <- n_base + j
      lp <- lp + stats::dnorm(con_p[idx], mean = 0, sd = 10, log = TRUE)
    }
  }
  idx_f <- (n_base + n_cov + 1):(n_base + n_cov + n_frail)
  if (frailty == "posstable") {
    lp <- lp + stats::dunif(con_p[idx_f[1]], 0, 1, log = TRUE)
  } else if (frailty == "mixture_invgauss") {
    lp <- lp + stats::dunif(con_p[idx_f[1]], 0, 1, log = TRUE) +
          stats::dgamma(con_p[idx_f[2]], shape = 0.01, rate = 0.01, log = TRUE) +
          stats::dgamma(con_p[idx_f[3]], shape = 0.01, rate = 0.01, log = TRUE)
  } else {
    for (k in seq_along(idx_f)) {
      lp <- lp + stats::dgamma(con_p[idx_f[k]], shape = 0.01, rate = 0.01, log = TRUE)
    }
  }
  lp
}
