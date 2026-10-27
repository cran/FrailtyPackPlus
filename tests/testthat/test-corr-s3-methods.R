test_that("S3 generic methods operate properly", {
  data(kidney_infection)
  expect_equal(nrow(kidney_infection), 38)
  
  fit <- corrFrailty(Surv(t1, d1, t2, d2) ~ age + sex, data = kidney_infection,
                     frailty = "invgaussian", baseline = "exponential",
                     method = "mcmc",
                     control = corrFrailty.control(n.iter = 50, n.burnin = 10,
                                                   n.thin = 1, n.chains = 1, seed = 123))
  
  # print and summary
  expect_output(print(fit))
  s <- summary(fit)
  expect_s3_class(s, "summary.corrFrailty")
  expect_output(print(s))
  
  # extractors
  cf <- coef(fit)
  expect_true(length(cf) > 0)
  v <- vcov(fit)
  expect_equal(dim(v), c(length(cf), length(cf)))
  ll <- logLik(fit)
  expect_s3_class(ll, "logLik")
  expect_true(is.numeric(AIC(fit)))
  expect_true(is.numeric(BIC(fit)))
  
  # predict & residuals
  pred_surv <- predict(fit, times = c(10, 50, 100), type = "survival")
  expect_equal(nrow(pred_surv), 3)
  expect_true(all(pred_surv$survival1 >= 0 & pred_surv$survival1 <= 1))
  
  res_cs <- residuals(fit, type = "cox-snell")
  expect_equal(nrow(res_cs), 38)
  
  # plot call
  expect_silent({
    pdf(NULL)
    plot(fit, which = c("survival", "hazard"))
    dev.off()
  })
})
