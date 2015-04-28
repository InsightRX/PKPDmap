library(Rcpp)
library(PKPDsim)
library(devtools)
install_github("ronkeizer/PKPDsim")
library(bbmle)
library(mvtnorm)

## define parameters
pk1 <- new_ode_model(model = "pk_1cmt_oral")
regimen  <- new_regimen (amt = 100, interval = 12, n = 2)
parameters   <- list("CL" = 5, "V" = 50, "KA" = 0.5) 
omega = cv_to_omega(list("CL" = 0.2, "V" = 0.2, "KA" = 0.1), p)

## simulate single individual in population
data <- sim_ode(ode = "pk1", parameters = list(CL = 6, V = 55, KA=0.1), t_obs = c(1, 4, 12, 18, 23))

## Fit
fit <- get_map_estimates(model = pk1, 
                         data = data, 
                         parameters = parameters,
                         omega = omega,
                         regimen = regimen,
                         error = list(prop = 0.1, add = 0.1))
fit

print.map_estimates <- function(obj) {
  for (i in seq(obj$parameters)) {
    cat(paste0(names(obj$parameters[i]), ": ", round(as.numeric(obj$parameters[i]), 3), "\n"))
  }
}