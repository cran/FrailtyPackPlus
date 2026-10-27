test_that("generate_powerset works as expected", {
  ps1 <- generate_powerset(1)
  expect_equal(length(ps1), 2)
  expect_equal(ps1[[1]], integer(0))
  expect_equal(ps1[[2]], 1)
  
  ps2 <- generate_powerset(2)
  expect_equal(length(ps2), 4)
  expect_equal(ps2[[1]], integer(0))
  expect_equal(ps2[[4]], c(1, 2))
  
  ps3 <- generate_powerset(3)
  expect_equal(length(ps3), 8)
  
  expect_error(generate_powerset(0))
  expect_error(generate_powerset(-1))
  expect_error(generate_powerset(NA))
})

test_that("compute_probabilities and loglik_oneshot work accurately", {
  data(classH, package = "FrailtyPackPlus")
  X_cov <- cbind(1, classH$stress)
  tau <- classH$time
  count_mat <- as.matrix(classH[, c("n_none", "n_T", "n_G", "n_both")])
  par_paper <- c(-4.9897, -0.0058, -16.4365, 0.0183, 0.5)
  
  prob_res <- compute_probabilities(par_paper, X_cov, tau, M = 2, baseline = "exponential", frailty = TRUE)
  expect_equal(nrow(prob_res$prob_matrix), 8)
  expect_equal(ncol(prob_res$prob_matrix), 4)
  expect_true(all(prob_res$prob_matrix >= 0 & prob_res$prob_matrix <= 1))
  expect_equal(rowSums(prob_res$prob_matrix), rep(1, 8), tolerance = 1e-6)
  
  ll <- loglik_oneshot(par_paper, count_mat, X_cov, tau, M = 2, baseline = "exponential", frailty = TRUE)
  expect_true(is.finite(ll))
  expect_true(ll < 0)
})
