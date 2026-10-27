test_that("sim_add_frailty generates valid datasets across all configurations", {
  d_gl <- sim_add_frailty(n = 20, beta = 0.15, frailty_type = "genlindley", baseline_dist = "gw", seed = 1)
  expect_equal(nrow(d_gl), 20)
  expect_true(all(d_gl$t1 > 0 & d_gl$t2 > 0))
  expect_true(all(d_gl$delta1 %in% c(0, 1) & d_gl$delta2 %in% c(0, 1)))

  d_wl <- sim_add_frailty(n = 20, beta = 0.15, frailty_type = "wlindley", baseline_dist = "gll2", seed = 2)
  expect_equal(nrow(d_wl), 20)
  expect_true(all(d_wl$t1 > 0 & d_wl$t2 > 0))

  d_gam <- sim_add_frailty(n = 20, beta = 0.15, frailty_type = "gamma", baseline_dist = "gw", seed = 3)
  expect_equal(nrow(d_gam), 20)

  d_ig <- sim_add_frailty(n = 20, beta = 0.15, frailty_type = "invgaussian", baseline_dist = "gll2", seed = 4)
  expect_equal(nrow(d_ig), 20)
})
