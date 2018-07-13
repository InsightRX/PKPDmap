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
omega <- c(0.0406, 
           0.0623, 0.117)
reg <- new_regimen(amt = 100000, times=c(0, 24), type="bolus")

fit <- get_map_estimates(parameters = par,
                         model = model,
                         regimen = new_regimen(amt = 100000, times=c(0, 24), type="bolus"),
                         omega = omega,
                         error = list(prop = 0, add = sqrt(1.73E+04)),
                         data = data1, 
                         mixture = list("CL" = c(5, 9)) )

assert("correct CL", round(fit$parameters$CL,3) == 8.235)
assert("correct V", round(fit$parameters$V,3) == 99.762)
assert("mixture object created", !is.null(fit$mixture$probabilities))
assert("mixture probs ok", all(round(fit$mixture$probabilities, 2) == c(0.82, 0.18)))

