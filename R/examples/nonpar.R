library(PKPDmap)
library(PKPDsim)
library(PKPDplot)
library(pk1cmtiv)
library(ggplot2)

###############################################################
## Initialize
###############################################################

model      <- pk1cmtiv::model()
parameters <- list(CL = 5, V = 50) # population
omega      <- c(0.1, 0.05, 0.1)
regimen    <- new_regimen(amt = 500, n = 3, interval = 12, type = "bolus")
ruv        <- list(prop = 0.1, add = 1)

t_obs <- c(24, 26, 36)
res <- sim(ode = model, parameters = parameters, regimen = regimen, 
           omega = omega, n = 100, only_obs = TRUE)

res %>% plot() + 
  geom_point(data = data.frame(t = t_obs, y = 0), (aes(x=t, y = y)))

###############################################################
## Simulate single outlying individual to be estimated
###############################################################

par1 <- list(CL = 10, V = 130)
res1 <- sim(ode = model, parameters = par1, regimen = regimen, only_obs = TRUE, t_obs = c(24, 26, 36))
res1$y <- res1$y # * (1 + rnorm(1, 0, ruv$prop)) + rnorm(1, 0, ruv$add)
fit <- get_map_estimates(data = data.frame(t = res1$t, y = res1$y, evid = 0), 
                         model = model, parameters = parameters, regimen = regimen, 
                         omega = omega, error = ruv)

fit_flat <- get_map_estimates(data = data.frame(t = res1$t, y = res1$y, evid = 0), 
                         model = model, parameters = parameters, regimen = regimen, 
                         omega = omega, error = ruv, weight_prior = 0.01)

fit_np1 <- get_map_estimates(data = data.frame(t = res1$t, y = res1$y, evid = 0), 
                            model = model, parameters = parameters, regimen = regimen, 
                            omega = omega, error = ruv,
                            type = "np_hybrid", np_settings = list(grid_adaptive = FALSE))

fit_np2 <- get_map_estimates(data = data.frame(t = res1$t, y = res1$y, evid = 0), 
                            model = model, parameters = parameters, regimen = regimen, 
                            omega = omega, error = ruv,
                            type = "np_hybrid", np_settings = list(grid_adaptive = TRUE))

###############################################################
## the fit will have some shrinkage (10-20%) as expected 
## now, we're going to put a grid around the MAP estimate
###############################################################

pars <- create_grid_around_parameters(fit$parameters, 
                                      span = .6, 
                                      exponential = TRUE,
                                      grid_size = 8) # get grid around MAP estimates

np <- get_np_estimates(parameter_grid = pars,
                         error = list(prop = ruv$prop/2, add=ruv$add/2), 
                         model = model,
                         regimen = regimen,
                         t_obs = t_obs,
                         data = res1$y)
np$parameters
