test_that("oneshot_frailty model fitting works on Class-H data", {
  data(classH, package = "FrailtyPackPlus")
  start_paper <- c(-3.997, -0.009, -14.412, 0.003, 0.004, -0.001, -5.493e-4, -0.005, 0.082)
  fit <- oneshot_frailty(
    data = classH,
    time = "time",
    stress = "stress",
    counts = c("n_none", "n_T", "n_G", "n_both"),
    components = c("Turn", "Ground"),
    baseline = "weibull",
    frailty = TRUE,
    start = start_paper,
    method = "mle"
  )
  
  expect_s3_class(fit, "oneshot_frailty")
  expect_equal(fit$M, 2)
  expect_equal(fit$baseline, "weibull")
  expect_true(fit$frailty)
  expect_equal(fit$loglik, -49.358, tolerance = 0.1)
  expect_equal(length(fit$coefficients), 9)
  expect_true(fit$coefficients["beta"] > 0)
  
  expect_output(print(fit))
  s <- summary(fit)
  expect_s3_class(s, "summary.oneshot_frailty")
  expect_output(print(s))
  
  cf <- coef(fit)
  expect_equal(length(cf), 9)
  
  vc <- vcov(fit)
  expect_equal(dim(vc), c(9, 9))
  
  ci <- confint(fit)
  expect_equal(dim(ci), c(9, 2))
  expect_true(ci["beta", "LCI"] >= 0)
  
  ll_obj <- logLik(fit)
  expect_equal(as.numeric(ll_obj), fit$loglik)
  
  expect_equal(AIC(fit), fit$aic)
  expect_equal(BIC(fit), fit$bic)
  
  fp <- fitted(fit, type = "probs")
  expect_equal(dim(fp), c(8, 4))
  
  fc <- fitted(fit, type = "counts")
  expect_equal(dim(fc), c(8, 4))
  
  res <- residuals(fit, type = "pearson")
  expect_equal(dim(res), c(8, 4))
})

test_that("oneshot_frailty works with exponential baseline and independence", {
  data(classH, package = "FrailtyPackPlus")
  fit_exp_ind <- oneshot_frailty(
    data = classH,
    time = "time",
    stress = "stress",
    counts = c("n_none", "n_T", "n_G", "n_both"),
    components = c("Turn", "Ground"),
    baseline = "exponential",
    frailty = FALSE,
    method = "mle"
  )
  
  expect_s3_class(fit_exp_ind, "oneshot_frailty")
  expect_equal(fit_exp_ind$baseline, "exponential")
  expect_false(fit_exp_ind$frailty)
  expect_equal(length(fit_exp_ind$coefficients), 4)
})

test_that("oneshot_frailty formula interface works", {
  data(classH, package = "FrailtyPackPlus")
  start_paper <- c(-3.997, -0.009, -14.412, 0.003, 0.004, -0.001, -5.493e-4, -0.005, 0.082)
  fit_form <- oneshot_frailty(
    formula = cbind(n_none, n_T, n_G, n_both) ~ stress,
    data = classH,
    time = "time",
    components = c("Turn", "Ground"),
    baseline = "weibull",
    frailty = TRUE,
    start = start_paper
  )
  expect_s3_class(fit_form, "oneshot_frailty")
  expect_equal(fit_form$n_params, 9)
})
