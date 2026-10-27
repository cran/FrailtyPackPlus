test_that("baseline hazard and cumulative hazard functions work accurately", {
  dists <- c("genweibull", "genloglogistic", "weibull", "exponential", "pareto", "lfr", "burr")
  t_vals <- c(0.5, 1.0, 2.5, 10.0)
  
  for (d in dists) {
    pnames <- corrFrailty:::baseline_param_names(d)
    params <- rep(0.5, length(pnames))
    if ("gamma" %in% pnames) params[which(pnames == "gamma")] <- 1.2
    
    res <- corrFrailty:::baseline_H_h(t_vals, d, params)
    expect_true(all(res$H > 0))
    expect_true(all(res$h > 0))
    expect_false(any(is.na(res$H)))
    expect_false(any(is.na(res$h)))
    
    # Inversion test: H_inv(H(t)) should equal t
    t_recov <- corrFrailty:::baseline_inv_H(res$H, d, params)
    expect_equal(t_recov, t_vals, tolerance = 1e-4)
  }
})
