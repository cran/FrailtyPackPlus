test_that("Generalized Weibull baseline functions work correctly", {
  t_vals <- c(5, 10, 20)
  delta <- 0.05
  zeta <- 2.0
  xi <- 0.7

  s <- s_gw(t_vals, delta, zeta, xi)
  ch <- cumh_gw(t_vals, delta, zeta, xi)
  h <- h_gw(t_vals, delta, zeta, xi)

  expect_true(all(s > 0 & s <= 1))
  expect_true(all(diff(s) < 0)) # Monotonically decreasing survival
  expect_equal(ch, -log(s), tolerance = 1e-6)
  expect_true(all(h > 0))
})

test_that("Generalized Log-Logistic baseline functions work correctly", {
  t_vals <- c(5, 10, 20)
  delta <- 0.05
  zeta <- 3.0
  xi <- 1.2

  s <- s_gll2(t_vals, delta, zeta, xi)
  ch <- cumh_gll2(t_vals, delta, zeta, xi)
  h <- h_gll2(t_vals, delta, zeta, xi)

  expect_true(all(s > 0 & s <= 1))
  expect_true(all(diff(s) < 0))
  expect_equal(ch, -log(s), tolerance = 1e-6)
  expect_true(all(h > 0))
})

test_that("Generalized Lindley Laplace transform and derivatives match numerical approximations", {
  s_val <- 0.5
  theta <- 2.5
  mu <- 1.5

  l0 <- laplace_genlindley(s_val, theta, mu)
  l1 <- laplace_deriv1_genlindley(s_val, theta, mu)
  l2 <- laplace_deriv2_genlindley(s_val, theta, mu)

  expect_true(l0 > 0 && l0 <= 1)
  expect_true(l1 < 0)
  expect_true(l2 > 0)

  # Numerical derivative check
  eps <- 1e-6
  l0_plus <- laplace_genlindley(s_val + eps, theta, mu)
  l0_minus <- laplace_genlindley(s_val - eps, theta, mu)
  l1_num <- (l0_plus - l0_minus) / (2 * eps)
  expect_equal(l1, l1_num, tolerance = 1e-5)

  l1_plus <- laplace_deriv1_genlindley(s_val + eps, theta, mu)
  l1_minus <- laplace_deriv1_genlindley(s_val - eps, theta, mu)
  l2_num <- (l1_plus - l1_minus) / (2 * eps)
  expect_equal(l2, l2_num, tolerance = 1e-5)

  v <- var_genlindley(theta, mu)
  expect_true(is.numeric(v) && v > 0)
})

test_that("Weighted Lindley Laplace transform and derivatives match numerical approximations", {
  s_val <- 0.5
  s_par <- 1.5

  l0 <- laplace_wlindley(s_val, s_par)
  l1 <- laplace_deriv1_wlindley(s_val, s_par)
  l2 <- laplace_deriv2_wlindley(s_val, s_par)

  expect_true(l0 > 0 && l0 <= 1)
  expect_true(l1 < 0)
  expect_true(l2 > 0)

  eps <- 1e-6
  l0_plus <- laplace_wlindley(s_val + eps, s_par)
  l0_minus <- laplace_wlindley(s_val - eps, s_par)
  l1_num <- (l0_plus - l0_minus) / (2 * eps)
  expect_equal(l1, l1_num, tolerance = 1e-5)

  v <- var_wlindley(s_par)
  expect_true(is.numeric(v) && v > 0)
})

test_that("Gamma and Inverse Gaussian Laplace transforms work correctly", {
  expect_equal(laplace_gamma(0, theta = 1.5), 1)
  expect_equal(laplace_invgaussian(0, theta = 1.5), 1)
  expect_equal(var_gamma(1.5), 1.5)
  expect_equal(var_invgaussian(1.5), 1.5)
})
