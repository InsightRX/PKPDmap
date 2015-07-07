library(Rcpp)
library(PKPDsim)
#library(devtools)
#install_github("ronkeizer/PKPDsim")
library(bbmle)
library(mvtnorm)
library(PKPDmap)

## define parameters
pk1 <- new_ode_model(model = "pk_1cmt_oral")
regimen  <- new_regimen (amt = 100, interval = 12, n = 2)
parameters   <- list("CL" = 5, "V" = 50, "KA" = 0.5) 
omega = cv_to_omega(list("CL" = 0.2, "V" = 0.2), parameters[1:2])

## simulate single individual in population
data <- sim_ode(ode = pk1, 
                parameters = list(CL = 6, V = 55, KA=0.4), 
                t_obs = c(0.5, 1, 1.5, 2, 4, 12, 18, 23),
                regimen = regimen)

## Fit
pk1 <- new_ode_model(model = "pk_1cmt_oral")
get_map_estimates(model = pk1, 
                  data = data, 
                  parameters = parameters,
#                  fixed = c("KA"),
                  omega = omega,
                  regimen = regimen,
                  error = list(prop = 0.1, add = 0.1))
