test_that("Bayesian MCMC and MLE fit successfully on simulated data", {
  set.seed(42)
  sim_dat <- sim_biv_rhr_frailty(n = 40, frailty = "WL", baseline = "GE", theta = 3.4, seed = 42)
  
  # MCMC
  fit_bayes <- fit_rhr_frailty_bayes(
    cbind(t1, t2, d1, d2) ~ age,
    data = sim_dat,
    frailty = "WL",
    baseline = "GE",
    n_iter = 50,
    burn_in = 10,
    thin = 1,
    verbose = FALSE
  )
  
  expect_s3_class(fit_bayes, "rhr_bayes_fit")
  expect_true(is.numeric(fit_bayes$estimates$SELF))
  expect_true(is.numeric(fit_bayes$estimates$ELF))
  expect_true(is.numeric(fit_bayes$estimates$PLF))
  expect_true(!is.na(fit_bayes$criteria$DIC))
  expect_true(!is.na(fit_bayes$kendall_tau))
  
  # MLE
  fit_mle <- fit_rhr_frailty_mle(
    cbind(t1, t2, d1, d2) ~ age,
    data = sim_dat,
    frailty = "WL",
    baseline = "GE"
  )
  
  expect_s3_class(fit_mle, "rhr_mle_fit")
  expect_true(is.numeric(fit_mle$coefficients))
  expect_true(!is.na(fit_mle$criteria$AIC))
  
  # Model comparison
  cmp <- compare_rhr_models(Bayes = fit_bayes, MLE = fit_mle)
  expect_equal(nrow(cmp), 2)
  expect_true(all(c("AIC", "BIC", "Kendall_Tau") %in% colnames(cmp)))
})
