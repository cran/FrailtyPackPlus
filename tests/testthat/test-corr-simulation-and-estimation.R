test_that("simulation and model estimation work seamlessly", {
  set.seed(42)
  sim_df <- sim_corrFrailty(n = 25, frailty = "gamma",
                            frailty.params = list(sigma = 0.5, rho = 0.5),
                            baseline = "weibull",
                            baseline.params = list(lambda = 0.05, gamma = 1.1),
                            beta = c(0.4), censor.rate = 0.15, seed = 42)
  
  expect_equal(nrow(sim_df), 25)
  expect_true(all(c("t1", "d1", "t2", "d2", "x1") %in% colnames(sim_df)))
  
  # MCMC fit
  fit_mcmc <- corrFrailty(Surv(t1, d1, t2, d2) ~ x1, data = sim_df,
                          frailty = "gamma", baseline = "weibull",
                          method = "mcmc",
                          control = corrFrailty.control(n.iter = 60, n.burnin = 10,
                                                        n.thin = 1, n.chains = 2, seed = 42))
  
  expect_s3_class(fit_mcmc, "corrFrailty")
  expect_equal(fit_mcmc$method, "mcmc")
  expect_true(is.numeric(fit_mcmc$coefficients))
  expect_false(any(is.na(fit_mcmc$coefficients)))
  expect_true(is.finite(fit_mcmc$logLik))
  expect_true(is.finite(fit_mcmc$DIC))
  expect_true(is.finite(fit_mcmc$LPML))
  
  # MLE fit
  fit_mle <- corrFrailty(Surv(t1, d1, t2, d2) ~ x1, data = sim_df,
                         frailty = "gamma", baseline = "weibull",
                         method = "mle",
                         control = corrFrailty.control(optim.control = list(maxit = 100)))
  
  expect_s3_class(fit_mle, "corrFrailty")
  expect_equal(fit_mle$method, "mle")
  expect_true(is.numeric(fit_mle$coefficients))
  expect_true(is.finite(fit_mle$logLik))
})
