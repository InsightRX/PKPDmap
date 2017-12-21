## Test MAP fit with LOQ data
library(testit)
library(PKPDmap)
library(PKPDsim)
Sys.setenv("R_TESTS" = "")
model <- new_ode_model("pk_1cmt_iv")
i <- 1

dat   <- read.table(file=paste0(system.file(package="PKPDmap"), "/pktab1"), skip=1, header=TRUE)
colnames(dat)[1:3] <- c("id", "t", "y")
data1 <- dat[dat$id == i & dat$EVID == 0,]
par   <- list(CL = 7.67, V = 97.7) 
lloq <- 500
data1$loq <- 0
data1[data1$y < lloq,]$loq <- -1
data1[data1$loq < 0,]$y  <- lloq
omega <- c(0.0406, 
           0.0623, 0.117)
reg <- new_regimen(amt = 100000, times=c(0, 24), type="bolus")

fit <- get_map_estimates(parameters = par,
                         model = model,
                         regimen = new_regimen(amt = 100000, times=c(0, 24), type="bolus"),
                         omega = omega,
                         error = list(prop = 0, add = sqrt(1.73E+04)),
                         data = data1,
                         censoring = "loq")
assert("correct CL", round(fit$parameters$CL,3) == 7.469)
assert("correct V", round(fit$parameters$V,3) == 92.382)

## Test ULOQ
uloq <- 1000
data1[data1$y > uloq,]$loq <- 1
data1[data1$loq > 0,]$y   <- uloq

fit <- get_map_estimates(parameters = par,
                         model = model,
                         regimen = new_regimen(amt = 100000, times=c(0, 24), type="bolus"),
                         omega = omega,
                         error = list(prop = 0, add = sqrt(1.73E+04)),
                         data = data1,
                         censoring = "loq")
assert("correct CL", round(fit$parameters$CL,2) == 7.54)
assert("correct V", round(fit$parameters$V,1) == 95.6)

