test_that("Loss function calculations SELF, MQSELF, PLF produce expected results", {
  set.seed(42)
  draws <- matrix(rgamma(500 * 2, shape = 5, rate = 2), ncol = 2)
  colnames(draws) <- c("p1", "p2")

  res <- calc_loss_estimates(draws)
  expect_equal(nrow(res), 2)
  expect_true(all(c("Est_SELF", "Risk_SELF", "Est_MQSELF", "Risk_MQSELF", "Est_PLF", "Risk_PLF") %in% colnames(res)))
  expect_true(all(res$Risk_SELF >= 0))
  expect_true(all(res$Risk_MQSELF >= 0))
  expect_true(all(res$Risk_PLF >= 0))
  expect_true(all(res$LCL <= res$UCL))
})
