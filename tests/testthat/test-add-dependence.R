test_that("Dependence measures Kendall's tau and frailty variances are valid", {
  kt_gl <- kendall_tau("genlindley", theta = 2.29, mu = 1.38)
  expect_true(kt_gl > 0 && kt_gl < 1)

  kt_wl <- kendall_tau("wlindley", s_param = 1.10)
  expect_true(kt_wl > 0 && kt_wl < 1)

  kt_gamma <- kendall_tau("gamma", theta = 1.5)
  expect_equal(kt_gamma, 1.5 / (1.5 + 2), tolerance = 1e-6)

  v_gl <- frailty_variance("genlindley", theta = 2.29, mu = 1.38)
  expect_true(v_gl > 0)

  v_wl <- frailty_variance("wlindley", s_param = 1.10)
  expect_true(v_wl > 0)
})
