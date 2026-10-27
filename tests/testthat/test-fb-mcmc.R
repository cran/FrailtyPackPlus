test_that("frailty_bayes executes adaptive MCMC and returns valid output", {
  data(kidney_frailty)
  fit <- frailty_bayes(Surv(time1, time2, status1, status2) ~ age,
                       data = kidney_frailty,
                       baseline = "weibull",
                       frailty = "gamma",
                       type = "multiplicative_bivariate",
                       n_iter = 100,
                       burn_in = 20,
                       thin = 1)
  
  expect_s3_class(fit, "frailty_bayes")
  expect_true(nrow(fit$chain) == 100)
  expect_true(is.finite(fit$diagnostics$DIC))
  expect_true(is.finite(fit$diagnostics$LPML))
  expect_true(fit$accept_rate > 0 && fit$accept_rate < 1)
  
  # Summary and print
  s <- summary(fit)
  expect_s3_class(s, "summary.frailty_bayes")
  expect_true(is.data.frame(fit$summary))
  expect_true(all(c("Mean", "SD", "HPD_Lower", "HPD_Upper") %in% names(fit$summary)))
})
