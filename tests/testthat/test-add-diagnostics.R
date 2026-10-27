test_that("MCMC convergence diagnostics and KS tests operate correctly", {
  set.seed(123)
  ch1 <- matrix(rnorm(500 * 2, mean = 2, sd = 0.2), ncol = 2)
  ch2 <- matrix(rnorm(500 * 2, mean = 2, sd = 0.2), ncol = 2)
  colnames(ch1) <- colnames(ch2) <- c("zeta1", "delta1")

  gew <- calc_geweke(ch1)
  expect_equal(nrow(gew), 2)
  expect_true(all(gew$p_value >= 0 & gew$p_value <= 1))

  gr <- calc_gelman_rubin(list(ch1, ch2))
  expect_equal(nrow(gr), 2)
  expect_true(all(gr$PSRF >= 0.9))

  hw <- calc_heidelberger_welch(ch1)
  expect_equal(nrow(hw), 2)
  expect_true(all(is.logical(hw$Stationary)))

  # KS test
  data(kidney_bivariate)
  ks_out <- calc_ks_test(kidney_bivariate$t1, kidney_bivariate$delta1, "gw", c(2.0, 0.05, 0.65))
  expect_true(ks_out$p_value >= 0 && ks_out$p_value <= 1)
})
