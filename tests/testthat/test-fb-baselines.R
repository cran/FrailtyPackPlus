test_that("eval_baseline works for all 10 baseline distributions", {
  t <- c(0.5, 1.0, 2.5, 10.0)
  
  # 1. weibull
  b_w <- eval_baseline(t, "weibull", c(0.1, 1.5))
  expect_equal(length(b_w$h0), length(t))
  expect_true(all(b_w$h0 > 0))
  expect_true(all(b_w$H0 > 0))
  expect_true(all(b_w$S0 >= 0 & b_w$S0 <= 1))
  
  # 2. exponential
  b_e <- eval_baseline(t, "exponential", c(0.2))
  expect_equal(b_e$h0, rep(0.2, length(t)))
  
  # 3. genexponential
  b_ge <- eval_baseline(t, "genexponential", c(0.1, 2.0))
  expect_true(all(b_ge$S0 >= 0 & b_ge$S0 <= 1))
  
  # 4. exppower
  b_ep <- eval_baseline(t, "exppower", c(0.1, 1.2))
  expect_true(all(b_ep$H0 >= 0))
  
  # 5. lfr
  b_lfr <- eval_baseline(t, "lfr", c(0.1, 0.05))
  expect_true(all(b_lfr$h0 > 0))
  
  # 6. genloglogistic
  b_gll <- eval_baseline(t, "genloglogistic", c(0.1, 1.2, 1.0))
  expect_true(all(b_gll$S0 >= 0 & b_gll$S0 <= 1))
  
  # 7. genweibull
  b_gw <- eval_baseline(t, "genweibull", c(0.01, 1.2, 0.5))
  expect_true(all(b_gw$S0 >= 0 & b_gw$S0 <= 1))
  
  # 8. genrayleigh
  b_gr <- eval_baseline(t, "genrayleigh", c(0.1, 1.5))
  expect_true(all(b_gr$S0 >= 0 & b_gr$S0 <= 1))
  
  # 9. weightedexp
  b_we <- eval_baseline(t, "weightedexp", c(0.1, 1.5))
  expect_true(all(b_we$S0 >= 0 & b_we$S0 <= 1))
  
  # 10. extendedweibull
  b_ew <- eval_baseline(t, "extendedweibull", c(0.01, 1.5))
  expect_true(all(b_ew$S0 >= 0 & b_ew$S0 <= 1))
})
