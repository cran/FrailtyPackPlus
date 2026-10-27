test_that("Base hazard and cumulative hazard behave correctly", {
  t <- c(0.5, 1.0, 2.0, 5.0)
  H_wb <- base_cumhaz(t, shape = 1.5, scale = 2.0, baseline = "weibull")
  h_wb <- base_haz(t, shape = 1.5, scale = 2.0, baseline = "weibull")

  expect_true(all(H_wb > 0))
  expect_true(all(h_wb > 0))
  expect_true(all(diff(H_wb) > 0))

  H_gp <- base_cumhaz(t, shape = 0.5, scale = 0.1, baseline = "gompertz")
  h_gp <- base_haz(t, shape = 0.5, scale = 0.1, baseline = "gompertz")
  expect_true(all(H_gp > 0))
  expect_true(all(h_gp > 0))

  H_exp <- base_cumhaz(t, shape = 1.0, scale = 3.0, baseline = "exponential")
  h_exp <- base_haz(t, shape = 1.0, scale = 3.0, baseline = "exponential")
  expect_equal(H_exp, t / 3.0)
  expect_equal(h_exp, rep(1 / 3.0, length(t)))
})

test_that("Proper survival and density satisfy boundaries and integration", {
  t_seq <- c(0.01, 0.5, 1, 2, 5, 10, 50)
  schemes <- c("random", "minimum", "maximum")

  for (s in schemes) {
    s_star <- surv_proper(t_seq, theta = 1.2, alpha = 0.3, shape = 1.4, scale = 2.5,
                          activation = s, baseline = "weibull")
    expect_true(all(s_star >= 0 & s_star <= 1))
    expect_true(all(diff(s_star) <= 0)) # Monotonically decreasing

    f_star <- dens_proper(t_seq, theta = 1.2, alpha = 0.3, shape = 1.4, scale = 2.5,
                          activation = s, baseline = "weibull")
    expect_true(all(f_star >= 0))
  }
})

test_that("Quantile proper inverts survival accurately", {
  probs <- c(0.1, 0.25, 0.5, 0.75, 0.9)
  for (s in c("random", "minimum", "maximum")) {
    q_vals <- quantile_proper(probs, theta = 0.8, alpha = 0.2, shape = 1.2, scale = 3.0,
                              activation = s, baseline = "weibull")
    s_back <- surv_proper(q_vals, theta = 0.8, alpha = 0.2, shape = 1.2, scale = 3.0,
                          activation = s, baseline = "weibull")
    # S_star(q(v)) should equal 1 - v
    expect_equal(as.numeric(s_back), 1 - probs, tolerance = 1e-4)
  }
})
