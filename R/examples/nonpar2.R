library(PKPDmap)
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

par1 <- list(CL = 4, V = 60)
res1 <- sim(ode = model, parameters = par1, 
            regimen = regimen, 
            only_obs = TRUE, 
            t_obs = c(24, 26, 36))
res1$y <- res1$y # * (1 + rnorm(1, 0, ruv$prop)) + rnorm(1, 0, ruv$add)

grid <- cbind(CL = rlnorm(50, log(5), .3), V = rlnorm(50, log(50), .3))
par <- get_np_estimates(parameter_grid = grid, 
                 error = ruv, model = model,
                 data = res1, regimen = regimen)
ggplot(par$prob, aes(x = CL, y = V, size = like)) + geom_point()

## Make plot of simulated lines based on estimated likelihoods
simdat <- c()
for (i in 1:length(par$prob[,1])) {
  tmp <- sim(ode = model, parameters = par$prob[i,], 
              regimen = regimen, 
              only_obs = TRUE, t_obs = seq(0,36,.5))
  simdat <- rbind(simdat, cbind(i, tmp, par$prob[i,]$like))
}
names(simdat) <- c("id", "id2", "t", "comp", "y", "like")
ggplot(simdat, aes(x = t, y = y, alpha = like, group=id)) + geom_line() +
  geom_point(data=res1, aes(x = t, y = y), alpha = 1, colour="red", size = 2)


## Now with an extended grid
grid2 <- rbind(grid,
               create_grid_around_parameters(list(CL = 5, V = 50), 
                                             span = 1, 
                                             exponential = TRUE,
                                             grid_size = 8))
par2 <- get_np_estimates(parameter_grid = grid2, 
                        error = ruv, model = model,
                        data = res1, regimen = regimen)
par2_map <- get_map_estimates(parameters = list(CL = 5, V = 50),
                              omega =  c(0.1, 0.05, 0.1),
                         error = ruv, model = model,
                         data = res1, regimen = regimen)
ggplot(par2$prob, aes(x = CL, y = V, size = like, colour=like)) + geom_point()

## Make plot of simulated lines based on estimated likelihoods
simdat2 <- c()
for (i in 1:length(par2$prob[,1])) {
  tmp <- sim(ode = model, parameters = par2$prob[i,], 
             regimen = regimen, 
             only_obs = TRUE, t_obs = seq(0,36,.5))
  simdat2 <- rbind(simdat2, cbind(i, tmp, par2$prob[i,]$like))
}
simdat2_map <- sim(ode = model, parameters = par2_map$parameters, 
                            regimen = regimen, 
                            only_obs = TRUE, t_obs = seq(0,36,.5))

names(simdat2) <- c("id", "id2", "t", "comp", "y", "like")
simdat2[simdat2$like < 1e-3,]$like <- 1e-3
ggplot(simdat2, 
       aes(x = t, y = y, alpha = like, group=id)) + geom_line() +
  geom_line(data=simdat2_map, aes(x = t, y = y), alpha = 1, colour='lightgreen', size = 2) + 
  geom_point(data=res1, aes(x = t, y = y), alpha = 1, colour="red", size = 2) +
  theme_plain()

                 