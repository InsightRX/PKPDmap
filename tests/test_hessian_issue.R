## Sometimes the Hessian can't be inverted and the vcov matrix is then not positive definite.
## In that case get_map_estimates() should just return the original omega matrix,
## and throw a warning.
Sys.setenv("R_TESTS" = "")

## The below code creates a TDM scenario that will result in such a situation,
## this is with the vancomycin Carreno model and a single TDM.
par <- list(
  V = 25.76,
  SCLSlope = 0.036,
  K12 = 2.29,
  K21 = 1.44,
  SCLInter = 0.18,
  TDM_INIT = 0
)
omega <- c(0.205590,
           0.000000, 0.308640,
           0.000000, 0.000000, 1.116760,
           0.000000, 0.000000, 0.000000, 1.443300,
           0.000000, 0.000000, 0.000000, 0.000000, 0.057068)

## Set up model and regimen

## Commented out to avoid taking on a dependency for this one use of clinPK, but
## left in so we track where the egfr value was derived from
# egfr <- clinPK::calc_egfr(age = 65, weight = 156, height = 180, scr = 1.49, sex = "male", method = "cockcroft_gault_adjusted")$value
egfr <- 75.080589758495
covs <- list(CRCL = PKPDsim::new_covariate(value = egfr))
mod <- PKPDsim::new_ode_model(
  code = "
    CLi = SCLInter + SCLSlope * (CRCL*16.667) \
    Vi = V \
    Qi = K12 * Vi \
    V2i = Qi / K21 \
    dAdt[0] = -(CLi/V)*A[0] - K12*A[0] + K21*A[1] \
    dAdt[1] = K12*A[0] - K21*A[1] \
    dAdt[2] = A[0]/V
  ",
  pk_code = "",
  parameters = par,
  omega_matrix = omega,
  obs = list(cmt = 1, scale="V"),
  dose = list(cmt = 1, bioav = 1),
  declare_variables = c("CLi", "Vi", "Qi", "V2i"),
  covariates = covs)
reg <- PKPDsim::new_regimen(
  amt = c(1000, 1000),
  times = c(0, 19.75),
  t_inf = 1.5,
  type = "infusion")
t_tdm <- 19.75 + 13.75
data <- data.frame(t = t_tdm, y = 8.8)

testit::assert("throws warning",
  testit::has_warning(
    fit <- PKPDmap::get_map_estimates(model = mod,
                             parameters = par,
                             covariates = covs,
                             regimen = reg,
                             data = data,
                             omega = omega,
                             error = list(prop = 0.1, add = 0.1),
                             fixed = "TDM_INIT")
  )
)
testit::assert("vcov matrix == input omega matrix", all(fit$vcov == omega))
