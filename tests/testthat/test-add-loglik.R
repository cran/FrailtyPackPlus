test_that("Log-likelihood evaluates without error across all model combinations", {
  data(kidney_bivariate)
  X <- as.matrix(kidney_bivariate[, c("age", "sex")])

  # 1. GL + GW
  pars_gl_gw <- c(2.0, 0.05, 0.65, 2.7, 0.06, 0.67, 2.3, 1.4, -0.1, -9.0)
  ll1 <- loglik_add_frailty(pars_gl_gw, kidney_bivariate$t1, kidney_bivariate$delta1,
                           kidney_bivariate$t2, kidney_bivariate$delta2,
                           X, frailty_type = "genlindley", baseline_dist = "gw")
  expect_true(is.finite(ll1))

  # 2. GL + GLL2
  pars_gl_gll2 <- c(3.8, 0.05, 1.2, 4.5, 0.06, 0.9, 2.3, 1.4, -0.1, -9.0)
  ll2 <- loglik_add_frailty(pars_gl_gll2, kidney_bivariate$t1, kidney_bivariate$delta1,
                            kidney_bivariate$t2, kidney_bivariate$delta2,
                            X, frailty_type = "genlindley", baseline_dist = "gll2")
  expect_true(is.finite(ll2))

  # 3. WL + GW
  pars_wl_gw <- c(2.0, 0.05, 0.65, 2.7, 0.06, 0.67, 1.1, -0.1, -9.0)
  ll3 <- loglik_add_frailty(pars_wl_gw, kidney_bivariate$t1, kidney_bivariate$delta1,
                            kidney_bivariate$t2, kidney_bivariate$delta2,
                            X, frailty_type = "wlindley", baseline_dist = "gw")
  expect_true(is.finite(ll3))

  # 4. WL + GLL2
  pars_wl_gll2 <- c(3.8, 0.05, 1.2, 4.5, 0.06, 0.9, 1.1, -0.1, -9.0)
  ll4 <- loglik_add_frailty(pars_wl_gll2, kidney_bivariate$t1, kidney_bivariate$delta1,
                            kidney_bivariate$t2, kidney_bivariate$delta2,
                            X, frailty_type = "wlindley", baseline_dist = "gll2")
  expect_true(is.finite(ll4))
})
