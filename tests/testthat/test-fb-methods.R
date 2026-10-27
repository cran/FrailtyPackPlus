test_that("S3 methods for frailty models work as expected", {
  data(kidney_frailty)
  fit <- frailty_mle(Surv(time1, time2, status1, status2) ~ age + sex,
                     data = kidney_frailty, baseline = "weibull", frailty = "gamma",
                     type = "multiplicative_bivariate")
  
  expect_s3_class(fit, "frailty_mle")
  expect_true(is.numeric(coef(fit)))
  expect_true(is.matrix(vcov(fit)))
  expect_s3_class(logLik(fit), "logLik")
  expect_true(is.numeric(AIC(fit)))
  expect_true(is.numeric(BIC(fit)))
  
  res_cs <- residuals(fit, type = "cox-snell")
  expect_equal(nrow(res_cs), nrow(kidney_frailty))
  
  res_mg <- residuals(fit, type = "martingale")
  expect_equal(nrow(res_mg), nrow(kidney_frailty))
  
  pred_s <- predict(fit, type = "survival", times = c(10, 50, 100))
  expect_equal(length(pred_s), 3)
  expect_true(all(pred_s >= 0 & pred_s <= 1))
})
