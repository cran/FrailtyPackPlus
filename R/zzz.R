#' @import stats
#' @import graphics
#' @import grDevices
#' @import utils
#' @import methods
#' @import survival
#' @importFrom MASS mvrnorm ginv
#' @importFrom numDeriv grad hessian
NULL

## Declare global variables to satisfy R CMD check
utils::globalVariables(c(
  "time1", "time2", "status1", "status2", "time", "status",
  "id", "density", "quantile", "median", "sd", "var",
  "lines", "par", "hist", "acf", "abline", "legend",
  "axis", "box", "polygon", "dev.new", "gray", "heat.colors",
  "points", "title", "mtext", "layout", "grid", "rgb",
  "Group", "Device", "Comp", "Stress", "Inspection_Time",
  "n_none", "n_T", "n_G", "n_both", "n_all",
  "Turn", "Ground", "zeta1", "delta1", "zeta2", "delta2",
  "beta", "theta", "alpha", "lambda", "k", "t", "s",
  "classH", "mice_tumor", "sim_data4", "insulin_gdm",
  "kidney_infection", "kidney_frailty", "kidney_bivariate",
  "retinopathy_frailty", "austrian_twin"
))
