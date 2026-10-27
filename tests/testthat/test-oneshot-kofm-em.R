test_that("em_oneshot_kofm executes and converges", {
  data(classH, package = "FrailtyPackPlus")
  X_cov <- cbind(1, classH$stress)
  tau <- classH$time
  count_mat <- as.matrix(classH[, c("n_none", "n_T", "n_G", "n_both")])
  
  fit_em <- em_oneshot_kofm(
    count_mat = count_mat,
    X_cov = X_cov,
    tau = tau,
    M = 2,
    baseline = "exponential",
    frailty = TRUE,
    max_iter = 30,
    tol = 1e-4
  )
  
  expect_type(fit_em, "list")
  expect_true("par" %in% names(fit_em))
  expect_equal(length(fit_em$par), 5)
  expect_true(is.finite(fit_em$loglik))
})
