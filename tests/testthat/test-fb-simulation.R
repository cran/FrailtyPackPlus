test_that("sim_frailty creates valid bivariate and univariate survival datasets", {
  # Bivariate
  d_biv <- sim_frailty(n = 60, baseline = "weibull", frailty = "gamma",
                       type = "multiplicative_bivariate", beta = c(0.3, -0.4))
  expect_equal(nrow(d_biv), 60)
  expect_true(all(c("time1", "time2", "status1", "status2", "x1", "x2") %in% names(d_biv)))
  expect_true(all(d_biv$time1 > 0 & d_biv$time2 > 0))
  
  # Univariate
  d_uni <- sim_frailty(n = 50, baseline = "exponential", frailty = "invgauss",
                       type = "multiplicative_univariate", beta = c(0.2, -0.3))
  expect_equal(nrow(d_uni), 50)
  expect_true(all(c("time", "status", "x1", "x2") %in% names(d_uni)))
  expect_true(all(d_uni$time > 0))
})
