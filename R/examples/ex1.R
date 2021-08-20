library(PKPDsim)
library(devtools)
#install_github("")
library(mvtnorm)
library(PKPDmap)

## define parameters
pk1 <- new_ode_model(code = "dAdt[1] = -(CL/V)*A[1] + rate;\n", obs = list(scale="V", cmt=1))
regimen  <- new_regimen(amt = 100, interval = 12, times=c(0,24), type="infusion")
parameters   <- list("CL" = 3, "V" = 10) 
# omega = PKPDsim::cv_to_omega(list("CL" = 0.2, "V" = 0.2, "KA" = 0.1), parameters)
omega <- PKPDsim::cv_to_omega(list("CL" = 0.2, "V" = 0.2), parameters[1:2])

## simulate single individual in population
data <- sim_ode(ode = pk1, 
                parameters = list(CL = 6, V = 14), 
                regimen = regimen,
                int_step_size = 0.1,
                t_obs = c(23, 25, 27, 32, 36))

## Fit
get_map_estimates(model = pk1, 
                  data = data, 
                  parameters = parameters,
                  omega = omega,
                  regimen = regimen,
                  error = list(prop = 0.1, add = 0.1))
