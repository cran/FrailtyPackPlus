library(survival)
library(FrailtyPackPlus)

test_that("sim_ziprom and rzicr_promotion generate valid datasets", {
  set.seed(99)
  d_rnd <- sim_ziprom(n = 80, activation = "random")
  expect_true(is.data.frame(d_rnd))
  expect_equal(nrow(d_rnd), 80)
  expect_true(all(c("time", "status", "x") %in% names(d_rnd)))
  expect_true(all(d_rnd$time >= 0))
  expect_true(all(d_rnd$status %in% c(0, 1)))

  # Check presence of zeros and censored
  expect_true(any(d_rnd$time == 0))
  expect_true(any(d_rnd$status == 0))

  d_min <- rzicr_promotion(n = 50, activation = "minimum")
  expect_equal(nrow(d_min), 50)

  d_max <- rzicr_promotion(n = 50, activation = "maximum")
  expect_equal(nrow(d_max), 50)
})
