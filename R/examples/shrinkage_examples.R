## Try to implement shrinkage control
## Goal is to minimize shrinkage

library(PKPDmap)
library(dplyr)
library(PKPDsim)
library(PKPDplot)
library(pk1cmtiv)
library(ggplot2)
library(vpc)

###############################################################
## Initialize
###############################################################

model      <- pk1cmtiv::model()
parameters <- list(CL = 5, V = 50) # population
omega      <- c(0.04 0.01, 0.04)
regimen    <- new_regimen(amt = 1500, n = 3, interval = 12, type = "bolus")
ruv        <- list(prop = 0.2, add = 1)

t_obs <- c(24, 26, 36)
res <- sim(ode = model, parameters = parameters, regimen = regimen, 
           omega = omega, n = 100, t_obs = seq(0, 36, .25), only_obs = TRUE)

res %>% plot() + 
  geom_point(data = data.frame(t = t_obs, y = 0), (aes(x=t, y = y)))

###############################################################
## Simulate single outlying individual to be estimated
###############################################################
par1 <- list(CL = 1.5, V = 15)
res1 <- sim(ode = model, parameters = par1, 
            regimen = regimen, 
            only_obs = TRUE, 
            t_obs = c(24))

fit1 <- get_map_estimates(model = model, 
                  parameters = parameters, 
                  omega = omega,
                  regimen = regimen, 
                  data = res1,
                  error = ruv)
fit1sim <- sim(ode = model, 
               parameters = fit1$parameters, 
            regimen = regimen, 
            only_obs = TRUE, 
            t_obs = seq(0, 36, .25))

res %>% plot() + 
  geom_line(data = fit1sim, (aes(x=t, y = y)), size = 2, colour="green") + 
  geom_point(data = res1, aes(x=t, y =y), colour="red")

## We see misfit (due to shrinkage). What do we do?
## Approaches:
## 2. Get another sample a.s.a.p.
## 3. Reduce dependence on prior
## 4. Automated, using shrinkage-control
## 5. Use NP ?
## 6. Use NP w/ extended grid

## Approach 1:
par2 <- list(CL = 1.5, V = 15)
res2 <- sim(ode = model, parameters = par2, 
            regimen = regimen, 
            only_obs = TRUE, 
            t_obs = c(24, 26, 36))

fit2 <- get_map_estimates(model = model, 
                          parameters = parameters, 
                          omega = omega,
                          regimen = regimen, 
                          data = res2,
                          error = ruv)
fit2sim <- sim(ode = model, 
               parameters = fit2$parameters, 
               regimen = regimen, 
               only_obs = TRUE, 
               t_obs = seq(0, 36, .25))

res %>% plot() + 
  geom_line(data = fit2sim, (aes(x=t, y = y)), size = 2, colour="green") + 
  geom_point(data = res2, aes(x=t, y =y), colour="red")

## Approach 3:
fit3 <- get_map_estimates(weight_prior = 0.25, 
                          model = model, 
                          parameters = parameters, 
                          omega = omega,
                          regimen = regimen, 
                          data = res1,
                          error = ruv)
fit3sim <- sim(ode = model, 
               parameters = fit3$parameters, 
               regimen = regimen, 
               only_obs = TRUE, 
               t_obs = seq(0, 36, .25))

res %>% plot() + 
  geom_line(data = fit3sim, (aes(x=t, y = y)), size = 2, colour="green") + 
  geom_point(data = res1, aes(x=t, y =y), colour="red")

## Approach 4
## Automated, using shrinkage control
#' @param shrinkage_control automatically control individual shrinkage. Not used when `NULL` (default). Suggested value e.g. `0.05`` for allowing 5% shrinkage (averaged over all parameters).
#' 
fit4 <- map_shrinkage_control(shrinkage_control = 0.05, 
                      model = model, 
                      parameters = parameters, 
                      omega = omega,
                      regimen = regimen, 
                      data = res1 %>% mutate(y= 25),
                      error = ruv) 
fit4$shrinkage_control

fit4sim <- sim(ode = model, 
               parameters = fit4$parameters, 
               regimen = regimen, 
               only_obs = TRUE, 
               t_obs = seq(0, 36, .25))

res %>% plot() + 
  geom_line(data = fit4sim, (aes(x=t, y = y)), size = 2, colour="green") + 
  geom_point(data = res1, aes(x=t, y =y), colour="red")


