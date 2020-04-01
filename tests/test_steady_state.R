library(PKPDmap)
library(testit)
library(PKPDsim)
Sys.setenv("R_TESTS" = "")

## Set up model and regimen
model <- new_ode_model("pk_2cmt_oral")
par_true <- list(CL=4.5, V=45, Q=3.5, V2=125, KA=0.6)
par <- list(CL=4, V=50, Q=3, V2=150, KA=0.5)
interval <- 12
reg <- new_regimen(amt = 1000, n = 30*2, interval = interval, type = "oral")
t_tdm <- max(reg$dose_times) + c(1, 3, 6, 8, 12)
tdm <- sim(model, parameters = par_true, regimen = reg, t_obs = t_tdm, only_obs=TRUE)

## plot
# library(ggplot2)
# dat <- sim(model, parameters = par_true, regimen = reg, only_obs=TRUE, t_obs = c(0:500))
# dat %>% 
#   ggplot(aes(x = t, y = y)) +
#   geom_line()

# fit with full dosing history
tmp <- get_map_estimates(
  parameters = par,
  model = model,
  data = tdm,
  regimen = reg,
  fixed = c("Q", "V2", "KA"),
  omega = c(0.1, 0.05, 0.1),
  error = list(prop = 0.1, add = 0.1)
)

# fit with steady state assumption
reg_ss <- new_regimen(amt = 1000, n = 2, interval = 12, type = "oral")
reg_ss$dose_times <- reg_ss$dose_times + max(reg$dose_times)
tmp_ss <- get_map_estimates(
  parameters = par,
  model = model,
  data = tdm,
  regimen = reg_ss,
  fixed = c("Q", "V2", "KA"),
  omega = c(0.1, 0.05, 0.1),
  error = list(prop = 0.1, add = 0.1),
  steady_state_analytic = list(
    f = "2cmt_oral",
    map = list(CL = "CL", V = "V", Q = "Q", V2 = "V2", KA = "KA"),
    auc = TRUE
  )
)

## tests
delta <- function(x,ref) (x-ref)/ref
assert("CL equal", delta(tmp_ss$parameters$CL, tmp$parameters$CL) < 0.001)
assert("V equal", delta(tmp_ss$parameters$V, tmp$parameters$V) < 0.001)

#############################
## Model with covariates
#############################
library(PKPDmap)
library(testit)
library(PKPDsim)

## Set up model and regimen
# covariates <- list(WT = new_covariate(c(70, 80), times = c(0, 120)))
covariates <- list(WT = new_covariate(c(70, 120), times = c(0, 120)))
model2 <- new_ode_model(code = "
        CLi = CL * pow(WT/70.0, 0.75)
        dAdt[1] = -KA*A[1]
        dAdt[2] = KA*A[1] - (CLi/V)*A[2]
      ", obs = list(cmt = 2, scale = "V"),
  declare_variables = "CLi",
  dose = list(cmt = 1), 
  covariates = covariates)

par_true <- list(CL=2.5, V=45, KA=0.5)
par <- list(CL=2, V=30, KA=0.5)
interval <- 24
reg <- new_regimen(amt = 1000, n = 30, interval = interval, type = "oral")
t_tdm <- max(reg$dose_times) + c(1, 3, 6, 8, 12)
tdm <- sim(model2, parameters = par_true, 
           covariates = covariates, 
           regimen = reg, t_obs = t_tdm, only_obs=TRUE)

## plot
# library(ggplot2)
# dat <- sim(model2, parameters = par,
#            covariates = covariates, regimen = reg, only_obs=TRUE, t_obs = c(0:500))
# dat %>%
#   ggplot(aes(x = t, y = y)) +
#   geom_line()

# fit with full dosing history extrapolated to steady state (manually)
tmp2 <- get_map_estimates(
  parameters = par,
  model = model2,
  data = tdm,
  regimen = reg,
  covariates = covariates, 
  fixed = c("KA"),
  omega = c(0.1, 0.05, 0.1),
  error = list(prop = 0.1, add = 0.1)
)

# fit with full dosing history extrapolated to steady state,
# but now using steady_state=TRUE in PKPDsim.
reg_ss <- new_regimen(
  amt = 1000, n = 1, interval = 24, type = "oral", 
  ss = TRUE, n_ss = 30)
tdm2 <- tdm
tdm2$t <- tdm2$t - max(reg$dose_times)
tmp2a <- get_map_estimates(
  parameters = par,
  model = model2,
  data = tdm2,
  regimen = reg_ss,
  covariates = covariates, 
  fixed = c("KA"),
  omega = c(0.1, 0.05, 0.1),
  error = list(prop = 0.1, add = 0.1)
)
## tests
assert("CL equal", delta(tmp2a$parameters$CL, tmp2$parameters$CL) < 0.001)
assert("V equal", delta(tmp2a$parameters$V, tmp2$parameters$V) < 0.001)

# fit with linear steady state equation
reg_ss <- new_regimen(amt = 1000, n = 1, interval = 24, type = "oral")
reg_ss$dose_times <- reg_ss$dose_times + max(reg$dose_times)
tmp2b <- get_map_estimates(
  parameters = par,
  model = model2,
  data = tdm,
  regimen = reg_ss,
  fixed = c("KA"),
  covariates = covariates, 
  omega = c(0.1, 0.05, 0.1),
  error = list(prop = 0.1, add = 0.1),
  steady_state_analytic = list(
    linear = TRUE,
    f = "1cmt_oral",
    map = list(CL = "CLi", V = "V", KA = "KA"),
    auc = TRUE
  )
)

## tests
assert("CL equal", delta(tmp2b$parameters$CL, tmp2$parameters$CL) < 0.001)
assert("V equal", delta(tmp2b$parameters$V, tmp2$parameters$V) < 0.001)
