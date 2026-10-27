test_that("oneshot_kofm model fitting works on Class-H data", {
  data(classH, package = "FrailtyPackPlus")
  fit <- oneshot_kofm(
    data = classH,
    time = "time",
    stress = "stress",
    counts = c("n_none", "n_T", "n_G", "n_both"),
    components = c("Turn", "Ground"),
    baseline = "exponential",
    frailty = TRUE,
    method = "mle"
  )
  
  expect_s3_class(fit, "oneshot_kofm")
  expect_equal(fit$M, 2)
  expect_equal(fit$baseline, "exponential")
  expect_true(fit$frailty)
  expect_true(is.finite(fit$loglik))
  expect_equal(length(fit$coefficients), 5)
  expect_true(fit$coefficients["beta"] > 0)
  
  expect_output(print(fit))
  s <- summary(fit)
  expect_s3_class(s, "summary.oneshot_kofm")
  expect_output(print(s))
  
  cf <- coef(fit)
  expect_equal(length(cf), 5)
  
  vc <- vcov(fit)
  expect_equal(dim(vc), c(5, 5))
  
  ci <- confint(fit)
  expect_equal(dim(ci), c(5, 2))
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

test_that("oneshot_kofm works with independence model and formula interface", {
  data(classH, package = "FrailtyPackPlus")
  fit_ind <- oneshot_kofm(
    formula = cbind(n_none, n_T, n_G, n_both) ~ stress,
    data = classH,
    time = "time",
    components = c("Turn", "Ground"),
    baseline = "exponential",
    frailty = FALSE,
    method = "mle"
  )
  
  expect_s3_class(fit_ind, "oneshot_kofm")
  expect_equal(fit_ind$baseline, "exponential")
  expect_false(fit_ind$frailty)
  expect_equal(length(fit_ind$coefficients), 4)
})
