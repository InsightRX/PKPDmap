## Test MAP fit using LTBS
library(testit)
library(PKPDmap)
library(PKPDsim)
Sys.setenv("R_TESTS" = "")
model <- new_ode_model("pk_1cmt_iv")
i <- 1

dat   <- read.table(file=paste0(system.file(package="PKPDmap"), "/pktab1"), skip=1, header=TRUE)
colnames(dat)[1:3] <- c("id", "t", "y")
dat   <- dat[dat$id <= 20,] # not necessary to do all 100 id's, if bug will be clear from 20 ids too.
data1 <- dat[dat$id == i & dat$EVID == 0,]
par   <- list(CL = 7.67, V = 97.7) 
lloq <- 500
data1$lloq <- 0
data1[data1$y < lloq,]$lloq <- 1
omega <- c(0.0406, 
           0.0623, 0.117)
reg <- new_regimen(amt = 100000, times=c(0, 24), type="bolus")
fit1 <- get_map_estimates(parameters = par,
                         model = model,
                         regimen = new_regimen(amt = 100000, times=c(0, 24), type="bolus"),
                         omega = omega,
                         error = list(prop = 0.15),
                         data = data1,
                         ltbs = FALSE,
                         censoring = list(flag = "lloq", limit = lloq, type = "lower")
)
fit2 <- get_map_estimates(parameters = par,
                          model = model,
                          regimen = new_regimen(amt = 100000, times=c(0, 24), type="bolus"),
                          omega = omega,
                          error = list(add = 0.15),
                          data = data1,
                          ltbs = TRUE,
                          censoring = list(flag = "lloq", limit = lloq, type = "lower")
)
assert("correct CL", round(fit2$parameters$CL,3) == 7.217)
assert("correct V", round(fit2$parameters$V,3) == 95.112)
