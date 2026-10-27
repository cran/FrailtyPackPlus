test_that("Baseline distributions evaluate properly", {
  t_vals <- c(5, 10, 20)
  
  # GE
  d_ge <- d_rhr_baseline(t_vals, "GE", c(zeta = 12.4, delta = 0.11))
  p_ge <- p_rhr_baseline(t_vals, "GE", c(zeta = 12.4, delta = 0.11))
  r_ge <- rhr_baseline(t_vals, "GE", c(zeta = 12.4, delta = 0.11))
  R_ge <- crhr_baseline(t_vals, "GE", c(zeta = 12.4, delta = 0.11))
  q_ge <- q_rhr_baseline(p_ge, "GE", c(zeta = 12.4, delta = 0.11))
  
  expect_true(all(d_ge > 0))
  expect_true(all(p_ge > 0 & p_ge < 1))
  expect_true(all(r_ge > 0))
  expect_true(all(R_ge > 0))
  expect_equal(q_ge, t_vals, tolerance = 1e-4)
  
  # EG
  d_eg <- d_rhr_baseline(t_vals, "EG", c(zeta = 9.1, delta = 0.10))
  p_eg <- p_rhr_baseline(t_vals, "EG", c(zeta = 9.1, delta = 0.10))
  r_eg <- rhr_baseline(t_vals, "EG", c(zeta = 9.1, delta = 0.10))
  R_eg <- crhr_baseline(t_vals, "EG", c(zeta = 9.1, delta = 0.10))
  q_eg <- q_rhr_baseline(p_eg, "EG", c(zeta = 9.1, delta = 0.10))
  
  expect_true(all(d_eg > 0))
  expect_true(all(p_eg > 0 & p_eg < 1))
  expect_true(all(r_eg > 0))
  expect_true(all(R_eg > 0))
  expect_equal(q_eg, t_vals, tolerance = 1e-4)
})

test_that("Laplace transforms and derivatives match mathematical properties", {
  s_vals <- c(0.1, 0.5, 1.0, 2.0)
  
  # WL
  L_wl <- laplace_frailty(s_vals, "WL", theta = 3.4)
  dL_wl <- dlaplace_frailty(s_vals, "WL", theta = 3.4)
  d2L_wl <- d2laplace_frailty(s_vals, "WL", theta = 3.4)
  
  expect_true(all(L_wl > 0 & L_wl <= 1))
  expect_true(all(dL_wl < 0))
  expect_true(all(d2L_wl > 0))
  
  # GL
  L_gl <- laplace_frailty(s_vals, "GL", theta = 3.0, mu = 1.5)
  dL_gl <- dlaplace_frailty(s_vals, "GL", theta = 3.0, mu = 1.5)
  d2L_gl <- d2laplace_frailty(s_vals, "GL", theta = 3.0, mu = 1.5)
  
  expect_true(all(L_gl > 0 & L_gl <= 1))
  expect_true(all(dL_gl < 0))
  expect_true(all(d2L_gl > 0))
  
  # Kendall's tau is positive and in [0, 1]
  tau_wl <- kendall_tau_rhr("WL", theta = 3.4)
  expect_true(tau_wl >= 0 && tau_wl <= 1)
  
  tau_gamma <- kendall_tau_rhr("Gamma", theta = 2.0)
  expect_equal(tau_gamma, 1 / 5, tolerance = 1e-4)
})
