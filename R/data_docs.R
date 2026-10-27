#' Dataset documentation for Frailty package


# --- From austrian_twin.R ---
#' Australian Twin Appendectomy Data
#'
#' Bivariate survival dataset containing ages at onset of appendectomy in 350 male twin pairs
#' from the Australian Twin Registry, with left-censoring.
#'
#' @format A data frame with 350 rows and 7 variables:
#' \describe{
#'   \item{t1}{Age at onset of appendectomy (in years) for twin 1 (lifetime).}
#'   \item{t2}{Age at onset of appendectomy (in years) for twin 2 (lifetime).}
#'   \item{d1}{Censoring indicator for twin 1 (1 = observed, 0 = left-censored).}
#'   \item{d2}{Censoring indicator for twin 2 (1 = observed, 0 = left-censored).}
#'   \item{age}{Age of the twin pair (covariate K1).}
#'   \item{k2}{Binary indicator of appendectomy condition for twin 1 (covariate K2).}
#'   \item{k3}{Binary indicator of appendectomy condition for twin 2 (covariate K3).}
#' }
#'
#' @source Duffy, D. L., Martin, N. G., & Mathews, J. D. (1990). Appendectomy in Australian twins.
#'   \emph{American Journal of Human Genetics}, 47(3), 590-592.
#'
#' @references
#' Tyagi, S., Pandey, A., & Chesneau, C. (2022). Identifying the effects of observed and unobserved
#' risk factors using weighted Lindley shared regression model. \emph{Journal of Statistical Theory and Practice},
#' 16(1), 16. \doi{10.1007/s42519-021-00241-9}.
#'
#' Tyagi, S., Pandey, A., & Chesneau, C. (2022). Weighted Lindley shared regression model for bivariate
#' left censored data. \emph{Sankhya B: The Indian Journal of Statistics}, 84(2), 754-784.
#' \doi{10.1007/s13571-022-00278-1}.
#'
#' Pandey, A., Hanagal, D. D., Gupta, P., & Tyagi, S. (2020). Analysis of Australian twin data using
#' generalized inverse Gaussian shared frailty models based on reversed hazard rate. \emph{International
#' Journal of Statistics and Reliability Engineering}, 7(2), 219-235.
#'
#' Pandey, A., Hanagal, D. D., Tyagi, S., & Gupta, P. (2021). Generalized Lindley shared frailty based
#' on reversed hazard rate. \emph{International Journal of Reliability, Quality and Safety Engineering},
#' 28(5), 2150040. \doi{10.1142/S0218539321500406}.
#'
#' Gupta, P., Pandey, A., & Tyagi, S. (2022). Comparison of multiplicative frailty models under generalized
#' log-logistic-II baseline distribution for counter heterogeneous left censored data. In \emph{Statistical
#' Techniques for Interdisciplinary Research} (pp. 97-108). ISBN: 978-93-91248-57-4.
#'
#' @examples
#' data(austrian_twin)
#' summary(austrian_twin)
"austrian_twin"



# --- From data.R ---
#' Modified Class-H Failure Modes Data
#'
#' @description
#' The modified Class-H failure modes data presented in Table 10 of Ling et al. (2021),
#' adapted from motorette insulation testing originally reported by Nelson (2009).
#' The dataset records the outcomes of one-shot inspections across 8 test groups
#' under four elevated temperature conditions with \eqn{M = 2} failure modes:
#' Turn failure (\code{T}) and Ground failure (\code{G}).
#'
#' @format A data frame with 8 rows and 7 variables:
#' \describe{
#'   \item{\code{group}}{Test group index (\eqn{i = 1, \dots, 8}).}
#'   \item{\code{stress}}{Accelerated temperature stress level in degrees Fahrenheit (\eqn{s_i = 374, 428, 464, 500}).}
#'   \item{\code{time}}{Pre-specified inspection time in hours (\eqn{\tau_i = 800, \dots, 10000}).}
#'   \item{\code{n_none}}{Number of functioning devices (neither component malfunctioned, \eqn{X = \emptyset}).}
#'   \item{\code{n_T}}{Number of devices with only Turn component failure (\eqn{X = \{T\}}).}
#'   \item{\code{n_G}}{Number of devices with only Ground component failure (\eqn{X = \{G\}}).}
#'   \item{\code{n_both}}{Number of devices with both Turn and Ground failures (\eqn{X = \{T, G\}}).}
#' }
#'
#' @source
#' Ling, M.-H., Balakrishnan, N., Yu, C., and So, H. Y. (2021). Inference for One-Shot
#' Devices with Dependent k-Out-of-M Structured Components under Gamma Frailty.
#' \emph{Mathematics}, 9(23), 3032. \doi{10.3390/math9233032}.
#'
#' Nelson, W. B. (2009). \emph{Accelerated Testing: Statistical Models, Test Plans,
#' and Data Analysis}. John Wiley & Sons, Hoboken, NJ.
#'
#' @examples
#' data(classH)
#' head(classH)
"classH"

#' ED01 Mice Tumor Toxicological Data
#'
#' @description
#' Toxicological study data from the ED01 experiment on 671 mice across 6 test groups
#' (Table 12 of Ling et al., 2021), originally discussed by Lindsey and Ryan (1993).
#' Animals were exposed to 2-acetylaminofluorene (2-AAF) at control (0 ppm) or high dose
#' (150 ppm) and inspected at 12, 18, or 33 months for \eqn{M = 2} endpoints:
#' Bladder Tumor onset (\code{T}) and Death (\code{D}).
#'
#' @format A data frame with 6 rows and 7 variables:
#' \describe{
#'   \item{\code{group}}{Test group index (\eqn{i = 1, \dots, 6}).}
#'   \item{\code{dose}}{Carcinogen 2-AAF dose level in ppm (\eqn{s_i = 0, 150}).}
#'   \item{\code{time}}{Sacrifice / inspection time in months (\eqn{\tau_i = 12, 18, 33}).}
#'   \item{\code{n_none}}{Number of sacrificed mice without tumor (\eqn{X = \emptyset}).}
#'   \item{\code{n_T}}{Number of sacrificed mice with bladder tumor (\eqn{X = \{T\}}).}
#'   \item{\code{n_D}}{Number of dead mice without bladder tumor (\eqn{X = \{D\}}).}
#'   \item{\code{n_both}}{Number of dead mice with bladder tumor (\eqn{X = \{D, T\}}).}
#' }
#'
#' @source
#' Ling, M.-H., Balakrishnan, N., Yu, C., and So, H. Y. (2021). Inference for One-Shot
#' Devices with Dependent k-Out-of-M Structured Components under Gamma Frailty.
#' \emph{Mathematics}, 9(23), 3032. \doi{10.3390/math9233032}.
#'
#' Lindsey, J. C. and Ryan, L. M. (1993). A three-state multiplicative model for rodent
#' tumorigenicity experiments. \emph{Journal of the Royal Statistical Society: Series C
#' (Applied Statistics)}, 42(2), 283--300. \doi{10.2307/2986234}.
#'
#' @examples
#' data(mice_tumor)
#' head(mice_tumor)
"mice_tumor"

#' Simulated One-Shot Device Test Data with 4 Failure Modes
#'
#' @description
#' A simulated one-shot device testing dataset with \eqn{M = 4} failure modes across
#' 6 test groups (Table 15 of Ling et al., 2021). Evaluated at three stress levels
#' (\eqn{s_i = 35, 45, 55}) and two inspection times (\eqn{\tau_i = 10, 20}) with
#' \eqn{K_i = 100} devices per group.
#'
#' @format A data frame with 6 rows and 19 variables:
#' \describe{
#'   \item{\code{group}}{Test group index (\eqn{i = 1, \dots, 6}).}
#'   \item{\code{stress}}{Accelerated stress level (\eqn{s_i = 35, 45, 55}).}
#'   \item{\code{time}}{Inspection time (\eqn{\tau_i = 10, 20}).}
#'   \item{\code{n_none}}{Number of devices with no malfunction (\eqn{X = \emptyset}).}
#'   \item{\code{n_1}}{Devices with only Component 1 failed.}
#'   \item{\code{n_2}}{Devices with only Component 2 failed.}
#'   \item{\code{n_3}}{Devices with only Component 3 failed.}
#'   \item{\code{n_4}}{Devices with only Component 4 failed.}
#'   \item{\code{n_12}}{Devices with Components 1 and 2 failed.}
#'   \item{\code{n_13}}{Devices with Components 1 and 3 failed.}
#'   \item{\code{n_14}}{Devices with Components 1 and 4 failed.}
#'   \item{\code{n_23}}{Devices with Components 2 and 3 failed.}
#'   \item{\code{n_24}}{Devices with Components 2 and 4 failed.}
#'   \item{\code{n_34}}{Devices with Components 3 and 4 failed.}
#'   \item{\code{n_123}}{Devices with Components 1, 2, and 3 failed.}
#'   \item{\code{n_124}}{Devices with Components 1, 2, and 4 failed.}
#'   \item{\code{n_134}}{Devices with Components 1, 3, and 4 failed.}
#'   \item{\code{n_234}}{Devices with Components 2, 3, and 4 failed.}
#'   \item{\code{n_all}}{Devices with all 4 components failed (\eqn{X = \Omega}).}
#' }
#'
#' @source
#' Ling, M.-H., Balakrishnan, N., Yu, C., and So, H. Y. (2021). Inference for One-Shot
#' Devices with Dependent k-Out-of-M Structured Components under Gamma Frailty.
#' \emph{Mathematics}, 9(23), 3032. \doi{10.3390/math9233032}.
#'
#' @examples
#' data(sim_data4)
#' head(sim_data4)
"sim_data4"



# --- From data_doc.R ---
#' Time to Insulin Initiation in Gestational Diabetes Mellitus Dataset
#'
#' A real-world dataset comprising retrospective times to insulin initiation for
#' 390 pregnant women diagnosed with gestational diabetes mellitus (GDM) before
#' 24 weeks of pregnancy, followed at the Hospital das Clinicas, University of Sao Paulo, Brazil.
#'
#' @format A data frame with 390 rows and 8 variables:
#' \describe{
#'   \item{\code{time}}{Time in weeks from the first obstetric visit to insulin initiation.
#'     Women initiating insulin within the first 7 days have time recorded as 0.}
#'   \item{\code{status}}{Event indicator: 1 = insulin initiated during pregnancy,
#'     0 = right-censored (cured/insulin-free throughout pregnancy).}
#'   \item{\code{family_diabetes}}{Family history of diabetes: 0 = No (n = 143), 1 = Yes (n = 247).}
#'   \item{\code{prior_gdm}}{Prior history of gestational diabetes: 0 = Primigravida (n = 94),
#'     1 = No prior GDM (n = 249), 2 = Yes (n = 47).}
#'   \item{\code{smoking}}{Smoking history: 0 = Non-smoker (n = 356), 1 = Smoker (n = 34).}
#'   \item{\code{prior_macrosomia}}{Prior fetal macrosomia: 0 = No (n = 356), 1 = Yes (n = 34).}
#'   \item{\code{hypertension}}{Prior chronic hypertension: 0 = No (n = 282), 1 = Yes (n = 108).}
#'   \item{\code{bmi_cat}}{Pre-pregnancy body mass index category:
#'     0 = Normal (BMI <= 25, n = 84), 1 = Overweight (25 < BMI <= 30, n = 120),
#'     2 = Obesity (BMI > 30, n = 186).}
#' }
#'
#' @source
#' Assuncao, D., Ramos, P. L., Tomazella, V., Leao, J., Osses, L., and Rodrigues, A. (2026).
#' Zero-Inflated Promotion Cure-Rate Models With Frailty and Flexible Activation Schemes.
#' \emph{Biometrical Journal}, 68:e70162. \doi{10.1002/bimj.70162}
#'
#' Calsavara, V. F., Rodrigues, A. S., Rocha, R., Louzada, F., et al. (2019).
#' Zero-Adjusted Defective Regression Models for Modeling Lifetime Data.
#' \emph{Journal of Applied Statistics}, 46(13), 2434-2459.
#'
#' @examples
#' data(insulin_gdm)
#' head(insulin_gdm)
#' table(insulin_gdm$family_diabetes)
#' table(insulin_gdm$prior_macrosomia)
#'
#' # Fit final selected model from paper with Random activation scheme
#' fit_random <- ziprom_fit(
#'   formula_p0 = ~ prior_macrosomia,
#'   formula_p1 = ~ smoking + family_diabetes,
#'   time = insulin_gdm$time,
#'   status = insulin_gdm$status,
#'   data = insulin_gdm,
#'   activation = "random",
#'   frailty = FALSE
#' )
#' summary(fit_random)
"insulin_gdm"



# --- From data.R ---
#' Kidney Infection Recurrence Bivariate Survival Data
#'
#' Recurrence times of infection at the catheter insertion point for 38 kidney
#' disease patients using portable dialysis machines, as studied in McGilchrist
#' and Aisbett (1991), Pandey et al. (2022), and Tyagi et al. (2024).
#'
#' @format A data frame with 38 rows (one per patient) and 10 variables:
#' \describe{
#'   \item{id}{Patient identification number (1 to 38).}
#'   \item{t1}{First infection recurrence time in days.}
#'   \item{delta1}{Censoring indicator for first recurrence (1 = infection occurred, 0 = censored).}
#'   \item{t2}{Second infection recurrence time in days.}
#'   \item{delta2}{Censoring indicator for second recurrence (1 = infection occurred, 0 = censored).}
#'   \item{age}{Patient age in years.}
#'   \item{sex}{Patient sex (0 = male, 1 = female).}
#'   \item{gn}{Indicator for Glomerulonephritis kidney disease (1 = yes, 0 = no).}
#'   \item{an}{Indicator for Acute Nephritis kidney disease (1 = yes, 0 = no).}
#'   \item{pkd}{Indicator for Polycystic Kidney Disease (1 = yes, 0 = no).}
#' }
#'
#' @source
#' McGilchrist, C. A., and Aisbett, C. W. (1991). Regression with frailty in
#' survival analysis. \emph{Biometrics}, 47(2), 461-466. \doi{10.2307/2532138}
#'
#' Pandey, A., Hanagal, D. D., and Tyagi, S. (2022). Generalised Lindley shared
#' additive frailty regression model for bivariate survival data.
#' \emph{Statistics in Transition new series}, 23(4), 161-176.
#' \doi{10.2478/stattrans-2022-0048}
#'
#' Tyagi, S., Pandey, A., Hanagal, D. D., and Chesneau, C. (2024). A study on
#' comparisons of additive regression frailty models to counter heterogeneity:
#' Bayesian strategies and case study. \emph{Communications in Statistics -
#' Simulation and Computation}. \doi{10.1080/03610918.2024.2391874}
#'
#' @examples
#' data(kidney_bivariate)
#' head(kidney_bivariate)
#' summary(kidney_bivariate)
"kidney_bivariate"



# --- From data.R ---
#' Kidney Infection Recurrence Data
#'
#' Recurrence times to kidney catheter infection for 38 kidney patients (McGilchrist and Aisbett, 1991).
#'
#' @format A data frame with 38 rows and 8 variables:
#' \describe{
#'   \item{id}{Patient identification number (1 to 38).}
#'   \item{time1}{Time to first catheter infection (days).}
#'   \item{time2}{Time to second catheter infection (days).}
#'   \item{status1}{Censoring indicator for first recurrence (1 = infection, 0 = censored).}
#'   \item{status2}{Censoring indicator for second recurrence (1 = infection, 0 = censored).}
#'   \item{age}{Age of patient in years.}
#'   \item{sex}{Gender of patient (1 = male, 0 = female).}
#'   \item{disease}{Disease type (GN, AN, PKD, Other).}
#' }
#' @source McGilchrist, C. A. and Aisbett, C. W. (1991). Regression with frailty in survival analysis. \emph{Biometrics}, 47: 461-466.
#' @references
#' Hanagal, D. D. and Dabade, A. D. (2013). A comparative study of shared frailty models for kidney infection data with generalized exponential baseline distribution. \emph{Journal of Data Science}, 11(1): 109-142.
"kidney_frailty"

#' Diabetic Retinopathy Study Paired Eye Data
#'
#' Bivariate survival times to vision loss in 197 diabetic patients with paired eyes.
#'
#' @format A data frame with 197 rows and 9 variables:
#' \describe{
#'   \item{id}{Patient identifier.}
#'   \item{time1}{Time to visual acuity loss for eye 1.}
#'   \item{time2}{Time to visual acuity loss for eye 2.}
#'   \item{status1}{Status indicator for eye 1 (1 = blindness, 0 = censored).}
#'   \item{status2}{Status indicator for eye 2 (1 = blindness, 0 = censored).}
#'   \item{laser}{Laser treatment type (1 = Xenon, 0 = Argon).}
#'   \item{age}{Age at diabetes diagnosis in years.}
#'   \item{risk}{Risk group indicator.}
#'   \item{trt}{Treatment status (0 = control eye, 1 = treated eye).}
#' }
#' @source Huster, W. J., Brookmeyer, R. and Self, S. G. (1989). Modelling paired survival data with covariates. \emph{Biometrics}, 45: 145-156.
"retinopathy_frailty"



# --- From data.R ---
#' Kidney Infection Bivariate Survival Dataset
#'
#' @description
#' Recurrence times to catheter infection for 38 kidney dialysis patients
#' from the classic study by McGilchrist and Aisbett (1991).
#' Each patient has two recorded recurrence times (first and second insertion of catheter).
#'
#' @format A \code{data.frame} with 38 rows (one per patient) and 9 variables:
#' \describe{
#'   \item{id}{Patient identification number (1 to 38).}
#'   \item{t1}{Time to first catheter infection or censoring (in days).}
#'   \item{d1}{Censoring indicator for first event (1 = infection, 0 = censored).}
#'   \item{t2}{Time to second catheter infection or censoring (in days).}
#'   \item{d2}{Censoring indicator for second event (1 = infection, 0 = censored).}
#'   \item{age}{Patient age in years.}
#'   \item{sex}{Patient sex (1 = male, 2 = female).}
#'   \item{GN}{Indicator for Glomerulonephritis disease type (1 = present, 0 = absent).}
#'   \item{AN}{Indicator for Acute Nephritis disease type (1 = present, 0 = absent).}
#'   \item{PKD}{Indicator for Polycystic Kidney Disease (1 = present, 0 = absent).}
#' }
#'
#' @source
#' McGilchrist, C. A. and Aisbett, C. W. (1991). Regression with frailty in survival analysis.
#' \emph{Biometrics}, 47(2), 461-466. <doi:10.2307/2532138>
#'
#' @references
#' Hanagal, D. D., Pandey, A., and Ganguly, A. (2017). Correlated gamma frailty
#' models for bivariate survival data. \emph{Communications in Statistics - Simulation
#' and Computation}, 46(5), 3627-3644. <doi:10.1080/03610918.2015.1085559>
#'
#' Hanagal, D. D. and Pandey, A. (2018). Correlated inverse Gaussian frailty
#' models for bivariate survival data. \emph{Communications in Statistics - Theory
#' and Methods}, 49(14), 3455-3475. <doi:10.1080/03610926.2018.1549256>
#'
#' @examples
#' data(kidney_infection)
#' head(kidney_infection)
"kidney_infection"


