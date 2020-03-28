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

## Compare steady state with full history:
# tdm_all <- sim(model, parameters = par_true, regimen = reg)
# tdm_all[tdm_all$t == (max(reg$dose_times) + interval),]
# f(1000, 12, parameters = par_true, covariates = NULL)

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
  steady_state = list(
    f = "TwoCompOral",
    map = list(CL = "CL", V = "V", Q = "Q", V2 = "V2", KA = "KA"),
    auc = TRUE
  )
)

## tests
delta <- function(x,ref) (x-ref)/ref
assert("CL equal", delta(tmp_ss$parameters$CL, tmp$parameters$CL) < 0.001)
assert("V equal", delta(tmp_ss$parameters$V, tmp$parameters$V) < 0.001)
