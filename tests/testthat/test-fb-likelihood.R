test_that("loglik_frailty evaluates correctly without NAs or NaNs", {
  data(kidney_frailty)
  
  # Multiplicative bivariate
  ll_mb <- frailty_mle(Surv(time1, time2, status1, status2) ~ age + sex,
                       data = kidney_frailty, baseline = "weibull", frailty = "gamma",
                       type = "multiplicative_bivariate", hessian = FALSE)
  expect_true(is.finite(ll_mb$logLik))
  expect_true(ll_mb$logLik < 0)
  
  # Additive bivariate
  ll_ab <- frailty_mle(Surv(time1, time2, status1, status2) ~ age + sex,
                       data = kidney_frailty, baseline = "weibull", frailty = "gamma",
                       type = "additive_bivariate", hessian = FALSE)
  expect_true(is.finite(ll_ab$logLik))
  expect_true(ll_ab$logLik < 0)
  
  # Univariate
  sim_u <- sim_frailty(n = 40, baseline = "weibull", frailty = "gamma",
                       type = "multiplicative_univariate", beta = 0.5)
  ll_u <- frailty_mle(Surv(time, status) ~ x1, data = sim_u,
                      baseline = "weibull", frailty = "gamma",
                      type = "multiplicative_univariate", hessian = FALSE)
  expect_true(is.finite(ll_u$logLik))
})
