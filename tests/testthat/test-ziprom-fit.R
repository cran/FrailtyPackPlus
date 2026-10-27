library(survival)
library(FrailtyPackPlus)

test_that("ziprom_fit fits correctly with simulated data across activation schemes", {
  set.seed(123)
  dat <- sim_ziprom(n = 100, beta1 = c(-1.5, 0.3), beta2 = c(-0.8, 0.5),
                    beta3 = c(0.1, 0), beta4 = c(1.0, 0), alpha = 0.25,
                    activation = "random", baseline = "weibull")

  # Random scheme with frailty
  fit_rnd <- ziprom_fit(Surv(time, status) ~ x, data = dat, activation = "random", frailty = TRUE)
  expect_s3_class(fit_rnd, "ziprom_fit")
  expect_true(is.numeric(fit_rnd$loglik))
  expect_true(!is.na(fit_rnd$aic))
  expect_true(all(fit_rnd$se > 0))
  expect_true(fit_rnd$p0_mean > 0 && fit_rnd$p0_mean < 1)
  expect_true(fit_rnd$p1_mean > 0 && fit_rnd$p1_mean < 1)

  # Minimum scheme without frailty
  fit_min <- ziprom_fit(Surv(time, status) ~ x, data = dat, activation = "minimum", frailty = FALSE)
  expect_s3_class(fit_min, "ziprom_fit")
  expect_equal(fit_min$frailty, FALSE)

  # Maximum scheme
  fit_max <- ziprom_fit(Surv(time, status) ~ x, data = dat, activation = "maximum", frailty = FALSE)
  expect_s3_class(fit_max, "ziprom_fit")

  # Gompertz baseline
  fit_gomp <- ziprom_fit(Surv(time, status) ~ x, data = dat, activation = "random",
                         baseline = "gompertz", frailty = FALSE)
  expect_s3_class(fit_gomp, "ziprom_fit")
})

test_that("S3 methods work as expected", {
  set.seed(42)
  dat <- sim_ziprom(n = 60, activation = "random")
  fit <- ziprom_fit(Surv(time, status) ~ x, data = dat, activation = "random")

  sm <- summary(fit)
  expect_s3_class(sm, "summary.ziprom_fit")
  expect_true(is.data.frame(sm$coef_table))

  cf <- coef(fit)
  expect_true(is.numeric(cf))
  expect_equal(length(cf), fit$npar)

  vc <- vcov(fit)
  expect_true(is.matrix(vc))
  expect_equal(dim(vc), c(fit$npar, fit$npar))

  ci <- confint(fit, level = 0.95)
  expect_true(is.matrix(ci))
  expect_equal(ncol(ci), 2)

  ll <- logLik(fit)
  expect_s3_class(ll, "logLik")
  expect_equal(AIC(fit), fit$aic)
  expect_equal(BIC(fit), fit$bic)

  r_q <- residuals(fit, type = "quantile")
  expect_equal(length(r_q), fit$n)

  r_cs <- residuals(fit, type = "cox-snell")
  expect_true(all(r_cs >= 0))
})
