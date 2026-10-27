library(survival)
library(FrailtyPackPlus)

test_that("Plotting methods run without error", {
  set.seed(42)
  dat <- sim_ziprom(n = 50, activation = "random")
  fit <- ziprom_fit(Surv(time, status) ~ x, data = dat, activation = "random")

  # Open a null graphics device
  pdf_null <- tempfile(fileext = ".pdf")
  grDevices::pdf(pdf_null)

  expect_silent(plot(fit, which = "survival"))
  expect_silent(plot(fit, which = "forecast"))
  expect_silent(plot(fit, which = "qq"))
  expect_silent(plot(fit, which = "residuals"))
  expect_silent(plot(fit, which = "hazard"))
  expect_silent(plot(fit, which = "density"))

  grDevices::dev.off()
  unlink(pdf_null)
})
