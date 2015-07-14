library(ggplot2)
library(PKPDsim)
#library(devtools)
#install_github("ronkeizer/PKPDsim")
library(bbmle)
library(mvtnorm)
library(PKPDmap)
library(pkbusulfanlt12kg)

## define parameters
pk1 <- pkbusulfanlt12kg::model("pk_busulfan_lt12kg")
regimen  <- new_regimen (amt = c(100,100,100,100), times=c(0, 24, 48, 72), t_inf=c(2,2,2,2), type="infusion")
parameters   <- list("CL" = 5, "V" = 50, 
                     "iov1" = 1, "iov2" = 1, "iov3" = 1, "iov4"= 1, 
                     "iov5" = 1, "iov6" = 1, "iov7" = 1, "iov8"= 1,
                     "iov9" = 1, "iov10" = 1) 
omega <- c(0.0622,
           0.0465, 0.0636,
           0, 0, 0.0169,
           0, 0, 0, 0.0169,
           0, 0, 0, 0, 0.0169, # rest of IOV is actaully dummy, not used.
           0, 0, 0, 0, 0, 0.0169,
           0, 0, 0, 0, 0, 0, 0.0169,
           0, 0, 0, 0, 0, 0, 0, 0.0169,
           0, 0, 0, 0, 0, 0, 0, 0, 0.0169,
           0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0169,
           0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0169,
           0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0169)

parameters   <- list("CL" = 5, "V" = 50, 
                     "iov1" = 1, "iov2" = 1) 
omega <- c(0.0622,
           0.0465, 0.0636,
           0, 0, 0.0169,
           0, 0, 0, 0.0169)

           
## simulate single individual in population
data <- sim_ode(ode = pk1, 
                omega = omega,
                parameters = parameters, 
                regimen = regimen)

ggplot(data[data$comp == "obs" & data$t < 72,], aes(x=t, y=y)) + geom_line()

## Fit
system.time({
  fit <- get_map_estimates(model = pk1, 
                    data = data, 
                    parameters = parameters,
                    method = "CG",
                    fixed = c("iov1", "iov2"),
                    # fixed = c("iov3", "iov4", "iov5", "iov6", "iov7", "iov8", "iov9", "iov10"),
                    omega = omega,
                    regimen = regimen,
                    int_step_size = 1,
                    error = list(prop = 0.1, add = 0.1))
})
