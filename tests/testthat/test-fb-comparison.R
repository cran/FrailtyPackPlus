test_that("model comparison functions operate correctly", {
  data(kidney_frailty)
  m1 <- frailty_bayes(Surv(time1, time2, status1, status2) ~ age,
                      data = kidney_frailty, frailty = "gamma",
                      n_iter = 100, burn_in = 20)
  m2 <- frailty_bayes(Surv(time1, time2, status1, status2) ~ age,
                      data = kidney_frailty, frailty = "invgauss",
                      n_iter = 100, burn_in = 20)
  
  # DIC, EAIC, EBIC, CPO, LPML
  d_val <- dic(m1)
  expect_true(is.finite(d_val$DIC))
  expect_true(is.finite(eaic(m1)))
  expect_true(is.finite(ebic(m1)))
  expect_true(is.finite(lpml(m1)))
  expect_equal(length(cpo(m1)$CPO), nrow(kidney_frailty))
  
  # Bayes Factor
  bf <- bayes_factor(m1, m2)
  expect_true(is.finite(bf$BF12))
  expect_true(is.character(bf$interpretation))
  
  # Compare table
  comp <- frailty_compare(m1, m2, names = c("Gamma", "InvGauss"))
  expect_equal(nrow(comp), 2)
  expect_true(all(c("DIC", "EAIC", "EBIC", "LPML") %in% names(comp)))
})
