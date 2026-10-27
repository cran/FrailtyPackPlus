test_that("fit_add_frailty_bayes and fit_add_frailty_mle estimate models properly", {
  data(kidney_bivariate)
  sub_data <- kidney_bivariate[1:15, ]

  # Bayesian fit
  fit_b <- fit_add_frailty_bayes(
    formula = ~ age,
    data = sub_data,
    frailty_type = "wlindley",
    baseline_dist = "gw",
    n_iter = 100,
    burn_in = 20,
    thin = 2,
    n_chains = 2,
    seed = 101
  )

  expect_s3_class(fit_b, "add_frailty_bayes")
  expect_s3_class(fit_b, "add_frailty")
  expect_true(is.data.frame(fit_b$estimates))
  expect_true(is.numeric(fit_b$kendall_tau))
  expect_true(is.numeric(fit_b$frailty_variance))

  # MLE fit
  fit_m <- fit_add_frailty_mle(
    formula = ~ age,
    data = sub_data,
    frailty_type = "wlindley",
    baseline_dist = "gw",
    control = list(maxit = 50)
  )

  expect_s3_class(fit_m, "add_frailty_mle")
  expect_s3_class(fit_m, "add_frailty")
  expect_true(is.numeric(coef(fit_m)))
  expect_true(is.matrix(vcov(fit_m)))
  expect_s3_class(logLik(fit_m), "logLik")
  expect_true(is.numeric(AIC(fit_m)))
  expect_true(is.numeric(BIC(fit_m)))

  # S3 methods
  expect_output(print(fit_b))
  s_b <- summary(fit_b)
  expect_output(print(s_b))

  expect_output(print(fit_m))
  s_m <- summary(fit_m)
  expect_output(print(s_m))

  pred_lp <- predict(fit_m, newdata = data.frame(age = c(30, 45)), type = "lp")
  expect_equal(length(pred_lp), 2)

  pred_s <- predict(fit_m, newdata = data.frame(age = 30), type = "survival", times = c(10, 20))
  expect_true(pred_s >= 0 && pred_s <= 1)
})
