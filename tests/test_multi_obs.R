## Example multiple observation types (with multiple residual errors)
library(testit)
library(PKPDsim)
library(PKPDmap)
Sys.setenv("R_TESTS" = "")

## define parameters
pk1 <- new_ode_model(code = "dAdt[1] = -(CL/V)*A[1]", obs = list(scale="V/1000", cmt=1))
regimen  <- new_regimen(amt = 100, interval = 12, n = 5, type="infusion", t_inf = 1)
parameters   <- list("CL" = 15, "V" = 150)
omega <- PKPDsim::cv_to_omega(list("CL" = 0.2, "V" = 0.2), parameters[1:2])

ruv_single <- list(prop = 0.1, add = 1)
ruv_multi <- list(prop = c(0.1, 1), add = c(0.1, 20))
ruv_multi2 <- list(prop = c(0.1, 1, 2), add = c(0.1, 1, 20))

## simulate single individual in population
# some observations with much higher residual error, should affect fit that much
set.seed(83475)
data_multi <- sim_ode(ode = pk1, 
                parameters = list(CL = 20, V = 200), 
                regimen = regimen,
                int_step_size = 0.1,
                only_obs = TRUE,
                obs_type =  c(1,2,1,2), 
                t_obs = c(2, 4, 6, 8),
                res_var = ruv_multi)

data_noruv <- sim_ode(ode = pk1, 
                      parameters = list(CL = 20, V = 200), 
                      regimen = regimen,
                      int_step_size = 0.1,
                      only_obs = TRUE,
                      t_obs = c(2, 4, 6, 8))

## Fit
# first fit with single low residual errror
fit1 <- get_map_estimates(model = pk1, 
                  data = data_multi, 
                  parameters = parameters,
                  omega = omega,
                  regimen = regimen,
                  obs_type_label = NULL,
                  error = ruv_single)

# then fit with multiple residual error. Fit should be affected less by points with large error
fit2 <- get_map_estimates(model = pk1, 
                  data = data_multi, 
                  parameters = parameters,
                  omega = omega,
                  regimen = regimen,
                  obs_type_label = "obs_type",
                  error = ruv_multi)
fit3 <- get_map_estimates(model = pk1, 
                          data = data_noruv, 
                          parameters = parameters,
                          omega = omega,
                          regimen = regimen,
                          error = ruv_single) 
testit::assert("Much more bias when high res.error points not weighted down", 
       (abs(fit3$parameters$CL - fit1$parameters$CL) / abs(fit3$parameters$CL - fit2$parameters$CL)) > 5)

## Test ordering of data
data_multi2 <- sim_ode(ode = pk1, 
                       parameters = list(CL = 20, V = 200), 
                       regimen = regimen,
                       int_step_size = 0.1,
                       only_obs = TRUE,
                       obs_type =  c(1, 2, 2, 1, 1, 2), 
                       t_obs = c(2, 4, 6, 6, 8, 8),
                       res_var = ruv_multi)
data_multi3 <- sim_ode(ode = pk1, 
                       parameters = list(CL = 20, V = 200), 
                       regimen = regimen,
                       int_step_size = 0.1,
                       only_obs = TRUE,
                       obs_type =  c(1, 2, 3, 2, 1, 1, 2), 
                       t_obs = c(2, 2, 2, 4, 4, 8, 8),
                       res_var = ruv_multi2)

# PKPDsim is expected to re-order. Testing here to make sure this is still the case
testit::assert(
  "PKPDsim re-orders as expected", 
  all(data_multi2$obs_type == c(1, 2, 1, 2, 1, 2))
)
testit::assert(
  "PKPDsim re-orders as expected", 
  all(data_multi3$obs_type == c(1, 2, 3, 1, 2, 1, 2))
)
fit4 <- get_map_estimates(model = pk1, 
                          data = data_multi3, 
                          parameters = parameters,
                          omega = omega,
                          regimen = regimen,
                          obs_type_label = "obs_type",
                          error = ruv_multi2)
testit::assert(
  "order is same for input vs output", 
  fit4$obs_type == data_multi3$obs_type
)
