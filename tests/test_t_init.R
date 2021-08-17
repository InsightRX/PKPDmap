## Test t_init functionality
## 
library(testit)
library(PKPDmap)
library(PKPDsim)
Sys.setenv("R_TESTS" = "")

## TDM before first dose:
## at t=-8, conc=10000
## Use this as true value in the simulations

dat   <- read.table(file=paste0(system.file(package="PKPDmap"), "/pktab1"), skip=1, header=TRUE)
colnames(dat)[1:3] <- c("id", "t", "y")
dat   <- dat[dat$id == 1 & dat$EVID == 0,]
par   <- list(CL = 7.67, V = 97.7, TDM_INIT = 500) 

mod <- new_ode_model(code = "
  dAdt[1] = -(CL/V)*A[1];
", state_init = "
  A[1] = TDM_INIT * V
", parameters = par, obs = list(cmt = 1, scale = "V"), cpp_show_code=F)
fits <- c()
omega <- c(0.0406, 
           0.0623, 0.117)
reg <- new_regimen(amt = 100000, times=c(0, 24), type="bolus")

fit4 <- get_map_estimates(parameters = par,
                         model = mod,
                         fixed = c("TDM_INIT"),
                         regimen = new_regimen(amt = 100000, times=c(0, 24), type="bolus"),
                         omega = omega,
                         weights = rep(1, length(dat$t)),
                         error = list(prop = 0, add = sqrt(1.73E+04)),
                         data = dat,
                         t_init = 4,
                         verbose = F
)

fit2 <- get_map_estimates(parameters = par,
                         model = mod,
                         fixed = c("TDM_INIT"),
                         regimen = new_regimen(amt = 100000, times=c(0, 24), type="bolus"),
                         omega = omega,
                         weights = rep(1, length(dat$t)),
                         error = list(prop = 0, add = sqrt(1.73E+04)),
                         data = dat,
                         t_init = 2,
                         verbose = F
)
fit0 <- get_map_estimates(parameters = par,
                         model = mod,
                         fixed = c("TDM_INIT"),
                         regimen = new_regimen(amt = 100000, times=c(0, 24), type="bolus"),
                         omega = omega,
                         weights = rep(1, length(dat$t)),
                         error = list(prop = 0, add = sqrt(1.73E+04)),
                         data = dat,
                         verbose = F
)
print(fit4)
print(str(fit4))
assert("check fit object returned", "map_estimates" %in% class(fit4))
assert("check fit object returned", "map_estimates" %in% class(fit2))
assert("check fit object returned", "map_estimates" %in% class(fit0))
assert(fit4$parameters$CL != fit2$parameters$CL)
assert(fit2$parameters$CL != fit0$parameters$CL)
