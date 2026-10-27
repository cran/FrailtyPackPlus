test_that("eval_laplace works for all 8 frailty distributions", {
  s <- c(0.0, 0.5, 1.2, 5.0)
  
  # 1. gamma
  fr_g <- eval_laplace(s, "gamma", 1.5)
  expect_equal(fr_g$L0[1], 1.0)
  expect_true(all(fr_g$L0 > 0 & fr_g$L0 <= 1))
  expect_true(all(fr_g$L1 > 0))
  expect_true(all(fr_g$L2 > 0))
  
  # 2. invgauss
  fr_ig <- eval_laplace(s, "invgauss", 1.0)
  expect_equal(fr_ig$L0[1], 1.0)
  expect_true(all(fr_ig$L0 > 0 & fr_ig$L0 <= 1))
  
  # 3. posstable
  fr_ps <- eval_laplace(s, "posstable", 0.7)
  expect_true(all(fr_ps$L0 > 0 & fr_ps$L0 <= 1))
  
  # 4. compoundpoisson
  fr_cp <- eval_laplace(s, "compoundpoisson", c(1.2, 1.0))
  expect_true(all(fr_cp$L0 > 0 & fr_cp$L0 <= 1))
  
  # 5. compoundnb
  fr_cnb <- eval_laplace(s, "compoundnb", c(1.0, 1.0, 1.0))
  expect_true(all(fr_cnb$L0 > 0 & fr_cnb$L0 <= 1))
  
  # 6. genlindley
  fr_gl <- eval_laplace(s, "genlindley", c(1.5, 0.5))
  expect_equal(fr_gl$L0[1], 1.0)
  expect_true(all(fr_gl$L0 > 0 & fr_gl$L0 <= 1))
  
  # 7. weightedlindley
  fr_wl <- eval_laplace(s, "weightedlindley", c(1.5, 1.2))
  expect_equal(fr_wl$L0[1], 1.0)
  expect_true(all(fr_wl$L0 > 0 & fr_wl$L0 <= 1))
  
  # 8. mixture_invgauss
  fr_mix <- eval_laplace(s, "mixture_invgauss", c(0.4, 1.0, 2.0))
  expect_equal(fr_mix$L0[1], 1.0)
  expect_true(all(fr_mix$L0 > 0 & fr_mix$L0 <= 1))
})

test_that("r_frailty generates positive draws", {
  n <- 50
  expect_equal(length(r_frailty(n, "gamma", 2.0)), n)
  expect_true(all(r_frailty(n, "gamma", 2.0) > 0))
  expect_true(all(r_frailty(n, "invgauss", 1.0) > 0))
  expect_true(all(r_frailty(n, "genlindley", c(1.0, 1.0)) > 0))
  expect_true(all(r_frailty(n, "weightedlindley", c(1.0, 1.5)) > 0))
})
