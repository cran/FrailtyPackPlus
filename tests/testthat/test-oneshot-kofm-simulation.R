test_that("sim_copula and sim_frailty generate valid failure times", {
  lam <- c(0.0086, 0.0067, 0.0052, 0.0024)
  
  sim_c <- sim_copula(n = 50, lambda = lam, beta = 0.3, seed = 123)
  expect_equal(dim(sim_c), c(50, 4))
  expect_true(all(sim_c >= 0))
  
  sim_f <- sim_frailty(n = 50, lambda = lam, beta = 0.3, seed = 123)
  expect_equal(dim(sim_f), c(50, 4))
  expect_true(all(sim_f >= 0))
})

test_that("sim_oneshot generates valid data frames", {
  a_par <- matrix(c(-6.0, 0.05, -6.5, 0.06), nrow = 2, byrow = TRUE)
  
  sim_d <- sim_oneshot(
    n_groups = 3,
    K = 40,
    inspection_times = c(10, 20, 30),
    stress_levels = c(35, 45, 55),
    a_params = a_par,
    beta = 0.2,
    seed = 123
  )
  
  expect_equal(nrow(sim_d), 3)
  count_cols <- c("n_none", "n_1", "n_2", "n_all")
  expect_true(all(count_cols %in% colnames(sim_d)))
  expect_equal(rowSums(sim_d[, count_cols]), rep(40, 3))
})
