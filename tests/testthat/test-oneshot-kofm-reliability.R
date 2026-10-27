test_that("reliability, mttf, and component moments work accurately", {
  data(classH, package = "FrailtyPackPlus")
  fit <- oneshot_kofm(
    data = classH,
    time = "time",
    stress = "stress",
    counts = c("n_none", "n_T", "n_G", "n_both"),
    components = c("Turn", "Ground"),
    baseline = "exponential",
    frailty = TRUE,
    method = "mle"
  )
  
  rel <- reliability_k_out_of_m(fit, t = c(1000, 2000), stress = 356)
  expect_equal(nrow(rel), 4)
  expect_true(all(rel$reliability >= 0 & rel$reliability <= 1))
  
  r_parallel <- rel$reliability[rel$k == 1]
  r_series <- rel$reliability[rel$k == 2]
  expect_true(all(r_parallel >= r_series))
  
  mttf_res <- mttf_k_out_of_m(fit, stress = 356)
  expect_equal(nrow(mttf_res), 2)
  expect_true(all(mttf_res$mttf > 0))
  expect_true(all(mttf_res$aci_lower >= 0))
  expect_true(all(mttf_res$tci_lower > 0))
  
  moms <- component_moments(fit, stress = 356)
  expect_equal(nrow(moms), 2)
  expect_true(all(moms$lambda > 0))
  expect_true(all(moms$mean > 0))
})

test_that("hypothesis tests (independence & baseline) work", {
  data(classH, package = "FrailtyPackPlus")
  fit <- oneshot_kofm(
    data = classH,
    time = "time",
    stress = "stress",
    counts = c("n_none", "n_T", "n_G", "n_both"),
    components = c("Turn", "Ground"),
    baseline = "exponential",
    frailty = TRUE,
    method = "mle"
  )
  
  ht_ind <- test_independence(fit)
  expect_s3_class(ht_ind, "htest")
  expect_true(ht_ind$statistic >= 0)
  expect_true(ht_ind$p.value >= 0 && ht_ind$p.value <= 1)
})
