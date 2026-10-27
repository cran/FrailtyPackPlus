library(survival)
library(FrailtyPackPlus)

test_that("predict and forecast_ziprom generate multi-level significance intervals", {
  set.seed(123)
  dat <- sim_ziprom(n = 70, activation = "random")
  fit <- ziprom_fit(Surv(time, status) ~ x, data = dat, activation = "random")

  # Predict survival
  t_eval <- c(0, 0.5, 1, 2, 5)
  s_pred <- predict(fit, type = "survival", times = t_eval)
  expect_true(all(s_pred >= 0 & s_pred <= 1))

  # Predict proper susceptible survival
  s_star_pred <- predict(fit, type = "susceptible_survival", times = t_eval)
  expect_true(all(s_star_pred >= 0 & s_star_pred <= 1))

  # Forecast table with 90%, 95%, 99% intervals
  fc <- forecast_ziprom(fit, times = t_eval)
  expect_true(is.data.frame(fc))
  expect_true(all(c("time", "survival", "se", "lower_90", "upper_90",
                    "lower_95", "upper_95", "lower_99", "upper_99") %in% names(fc)))

  # Check ordering of confidence bands: 99% wider than 95% wider than 90%
  expect_true(all(fc$lower_99 <= fc$lower_95 + 1e-6))
  expect_true(all(fc$lower_95 <= fc$lower_90 + 1e-6))
  expect_true(all(fc$upper_99 >= fc$upper_95 - 1e-6))
  expect_true(all(fc$upper_95 >= fc$upper_90 - 1e-6))

  # Evaluation metrics (MSE, RMSE, MAE, Brier Score, C-index)
  evals <- evaluate_ziprom(fit)
  expect_true(is.list(evals))
  expect_true(evals$MSE >= 0)
  expect_true(evals$RMSE >= 0)
  expect_true(evals$MAE >= 0)
  expect_true(evals$Brier_Score >= 0)
  expect_true(evals$C_index >= 0 && evals$C_index <= 1)
})
