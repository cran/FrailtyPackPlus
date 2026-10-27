test_that("frailty log-derivatives satisfy regularity conditions", {
  frailties <- list(
    list(type = "gamma", params = c(0.6, 0.4)),
    list(type = "invgaussian", params = c(0.6, 0.4)),
    list(type = "posstable", params = c(0.7, 0.4)),
    list(type = "comppoisson", params = c(0.5, 1.0, 1.0)),
    list(type = "none", params = numeric(0))
  )
  
  s1 <- c(0.1, 0.5, 1.5, 3.0)
  s2 <- c(0.2, 0.8, 1.2, 2.5)
  
  for (f in frailties) {
    pder <- corrFrailty:::log_phi_derivatives(s1, s2, f$type, f$params)
    
    expect_true(all(pder$phi <= 0), info = paste("phi <= 0 for", f$type))
    expect_true(all(pder$phi1 < 0), info = paste("phi1 < 0 for", f$type))
    expect_true(all(pder$phi2 < 0), info = paste("phi2 < 0 for", f$type))
    expect_true(all(pder$phi12 >= 0), info = paste("phi12 >= 0 for", f$type))
    expect_true(all(pder$phi12 + pder$phi1 * pder$phi2 > 0),
                info = paste("cross term > 0 for", f$type))
    expect_false(any(is.na(pder$phi)))
  }
})
