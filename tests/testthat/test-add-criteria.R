test_that("Information criteria and Bayes Factor work correctly", {
  crit <- calc_model_criteria(loglik = -300, k = 10, n = 50, deviance_samples = rnorm(100, 600, 5), deviance_at_mean = 600)
  expect_equal(crit$AIC, 2 * 10 - 2 * (-300))
  expect_equal(crit$BIC, 10 * log(50) - 2 * (-300))
  expect_true(is.finite(crit$AICc))
  expect_true(is.finite(crit$HQIC))
  expect_true(is.finite(crit$DIC))

  fit_dummy1 <- list(criteria = list(LogLik = -300), loglik = -300, k = 10, n = 50)
  fit_dummy2 <- list(criteria = list(LogLik = -305), loglik = -305, k = 10, n = 50)
  bf <- calc_bayes_factor(fit_dummy1, fit_dummy2)
  expect_true(is.numeric(bf) && bf > 0)
})
