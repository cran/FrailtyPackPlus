test_that("hypothesis tests (independence & baseline) work", {
  data(classH, package = "FrailtyPackPlus")
  start_paper <- c(-3.997, -0.009, -14.412, 0.003, 0.004, -0.001, -5.493e-4, -0.005, 0.082)
  fit_w <- oneshot_frailty(
    data = classH,
    time = "time",
    stress = "stress",
    counts = c("n_none", "n_T", "n_G", "n_both"),
    baseline = "weibull",
    frailty = TRUE,
    start = start_paper
  )
  
  ht_ind <- test_independence(fit_w)
  expect_s3_class(ht_ind, "htest")
  expect_true(ht_ind$statistic >= 0)
  expect_true(ht_ind$p.value >= 0 && ht_ind$p.value <= 1)
  
  ht_base <- test_baseline(fit_w)
  expect_s3_class(ht_base, "htest")
  expect_true(ht_base$statistic >= 0)
  expect_true(ht_base$p.value >= 0 && ht_base$p.value <= 1)
})

test_that("reliability, mttf, and component moments work", {
  data(classH, package = "FrailtyPackPlus")
  start_paper <- c(-3.997, -0.009, -14.412, 0.003, 0.004, -0.001, -5.493e-4, -0.005, 0.082)
  fit_w <- oneshot_frailty(
    data = classH,
    time = "time",
    stress = "stress",
    counts = c("n_none", "n_T", "n_G", "n_both"),
    baseline = "weibull",
    frailty = TRUE,
    start = start_paper
  )
  
  rel <- reliability_k_out_of_m(fit_w, t = c(1000, 2000), stress = 400)
  expect_equal(nrow(rel), 4)
  expect_true(all(rel$reliability >= 0 & rel$reliability <= 1))
  
  r_parallel <- rel$reliability[rel$k == 1]
  r_series <- rel$reliability[rel$k == 2]
  expect_true(all(r_parallel >= r_series))
  
  mttf_res <- mttf_k_out_of_m(fit_w, stress = 400, upper_limit = 50000)
  expect_equal(nrow(mttf_res), 2)
  expect_true(all(mttf_res$mttf > 0))
  
  moms <- component_moments(fit_w, stress = 400)
  expect_equal(nrow(moms), 2)
  expect_true(all(moms$lambda > 0))
})

test_that("sim_oneshot generates valid data frames", {
  a_par <- matrix(c(-6.0, 0.05, -6.3, 0.06), nrow = 2, byrow = TRUE)
  b_par <- matrix(c(-1.2, 0.03, -1.4, 0.03), nrow = 2, byrow = TRUE)
  
  sim_d <- sim_oneshot(
    n_groups = 3,
    K = 40,
    inspection_times = c(1000, 2000, 3000),
    stress_levels = c(30, 40, 50),
    a_params = a_par,
    b_params = b_par,
    beta = 0.1,
    baseline = "weibull",
    seed = 123
  )
  
  expect_equal(nrow(sim_d), 3)
  count_cols <- c("n_none", "n_1", "n_2", "n_both")
  expect_true(all(count_cols %in% colnames(sim_d)))
  expect_equal(rowSums(sim_d[, count_cols]), rep(40, 3))
})
