library(ggplot2)
library(dplyr)
library(PKPDsim)
library(devtools)
#install_github("")
library(bbmle)
library(mvtnorm)
library(PKPDmap)
library(pk1cmtiv)

## define parameters
# pk1 <- new_ode_model(code = "dAdt[1] = -(CL/V)*A[1]")
pk1 <- pk1cmtiv::model()
regimen  <- new_regimen(amt = 500, interval = 6, times=c(0,24), n = 5, type="infusion")
parameters   <- list("CL" = 5, "V" = 100) 
par_prior <- list("CL" = 7, "V" = 80)
# omega = PKPDsim::cv_to_omega(list("CL" = 0.2, "V" = 0.2, "KA" = 0.1), parameters)
omega <- PKPDsim::cv_to_omega(list("CL" = 0.4, "V" = 0.4), parameters[1:2])
ruv <- list(prop = 0.15, add = 0.2)

## simulate single individual in population
data <- sim_ode(ode = pk1, 
                parameters = parameters, 
                regimen = regimen,
                int_step_size = 0.1,
                only_obs = TRUE,
                res_var = ruv,
                t_obs = c(23, 25, 36, 45))

# # ## Fit
fit1 <- get_map_estimates(model = pk1,
                  data = data,
                  parameters = par_prior,
                  omega = omega,
                  regimen = regimen,
                  error = list(prop = 0.1, add = 0.5))
ipred1 <- sim_ode(ode = pk1,
                parameters = fit1$parameters,
                regimen = regimen,
                only_obs = TRUE,
                int_step_size = 0.1)
## Fit
ipred1 %>%
  filter(t<50) %>%
  ggplot(aes(x=t, y=y)) + geom_line() + geom_point(data=data) + scale_y_log10()

## Iterative fit
fits <- run_sequential_map(model = pk1, 
                          data = data, 
                          parameters = par_prior,
                          omega = omega,
                          regimen = regimen,
                          error = ruv)

fits$obs %>%
  filter(t < 50) %>%
  ggplot(aes(x = t, y = y)) +
    geom_ribbon(aes(ymin = y-1.96*stdv, ymax = y+1.96*stdv), fill="#bfbfbf") +
    geom_line() + 
    geom_line(data=ipred1, colour='darkblue', linetype = 'dashed') + 
    geom_point(data = data) +
    scale_y_log10()

fits$obs %>%
  filter(t %in% data$t) 
ipred1 %>%
  filter(t %in% data$t) 
