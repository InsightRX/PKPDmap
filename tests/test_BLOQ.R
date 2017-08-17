## Test MAP fit with LOQ data
library(testit)
library(PKPDmap)
library(PKPDsim)
Sys.setenv("R_TESTS" = "")

dat   <- read.table(file=paste0(system.file(package="PKPDmap"), "/pktab1"), skip=1, header=TRUE)
colnames(dat)[1:3] <- c("id", "t", "y")
dat   <- dat[dat$id <= 20,] # not necessary to do all 100 id's, if bug will be clear from 20 ids too.
par   <- list(CL = 7.67, V = 97.7) 
lloq <- 500
dat$lloq <- 0
dat[dat$y < lloq & dat$EVID == 0,]$lloq <- 1

model <- new_ode_model("pk_1cmt_iv")
fits <- c()
i <- 1
omega <- c(0.0406, 
           0.0623, 0.117)
reg <- new_regimen(amt = 100000, times=c(0, 24), type="bolus")
fit <- get_map_estimates(parameters = par,
                         model = model,
                         regimen = new_regimen(amt = 100000, times=c(0, 24), type="bolus"),
                         omega = omega,
                         weights = rep(1, length(dat[dat$id == i & dat$EVID == 0,1])),
                         error = list(prop = 0, add = sqrt(1.73E+04)),
                         data = dat[dat$id == i & dat$EVID == 0,],
                         censoring = list(flag = "lloq", limit = lloq, type = "lower") 
                         #residuals = TRUE, verbose = TRUE
)
