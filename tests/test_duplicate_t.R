## Test if datasets with duplicate measurements work.
## 
## Rationale:
## - users might want to use this to put more weight on samples (preferrable is to use weighting functionality, obviously)
## - data (e.g. from EMR) might have duplicates 
## We should show a warning that there are duplicates, but the estimation should not fail.

library(testit)
library(PKPDmap)
library(PKPDsim)
Sys.setenv("R_TESTS" = "")

dat   <- read.table(file=paste0(system.file(package="PKPDmap"), "/pktab1"), skip=1, header=TRUE)
colnames(dat)[1:3] <- c("id", "t", "y")
dat   <- dat[dat$id == 1 & dat$EVID == 0,]
par   <- list(CL = 7.67, V = 97.7) 

model <- new_ode_model("pk_1cmt_iv")
fits <- c()
omega <- c(0.0406, 
           0.0623, 0.117)
reg <- new_regimen(amt = 100000, times=c(0, 24), type="bolus")

## add duplicates at t=4 and t=24
## make it difficult, becasue we'll put t=24 exactly at same timepoint as new dose
##
dat[dat$t == 23.9,]$t <- 24 
dat <- rbind(dat, dat[dat$t %in% c(4, 24),])
dat <- dat[order(dat$t),]

fit <- get_map_estimates(parameters = par,
                         model = model,
                         regimen = new_regimen(amt = 100000, times=c(0, 24), type="bolus"),
                         omega = omega,
                         weights = rep(1, length(dat$t)),
                         error = list(prop = 0, add = sqrt(1.73E+04)),
                         data = dat
)
assert("check fit object returned", "map_estimates" %in% class(fit))
