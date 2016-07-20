library(testit)
library(PKPDmap)
library(PKPDsim)
library(dplyr)

Sys.setenv("R_TESTS" = "")

dat   <- read.table(file=paste0(system.file(package="PKPDmap"), "/pktab1"), skip=1, header=TRUE)
colnames(dat)[1:3] <- c("id", "t", "y")
# dat   <- dat[dat$id <= 20,] # not necessary to do all 100 id's, if bug will be clear from 20 ids too.
# ebe   <- ebe[ebe$id <= 20,]  
parameters   <- list(CL = 3, V = 40) # NONMEM estimates: list(CL = 7.67, V = 97.7) 
omega <- c(0.3,
           0.1, 0.3)
err <- 100
model <- new_ode_model("pk_1cmt_iv")
regimen <- new_regimen(amt = 100000, times=c(0, 24), type="bolus")

run <- FALSE

if(run) {
  run_its(parameters = parameters,
          omega = omega,
          err = 100,
          regimen = regimen,
          model = model,
          max_iter=10,
          data = dat[dat$id < 20,])
}

test <- c(0.01010668, 0.02881619, 0.08216081)

check_high_corr <- function(cov_mat, limit = 0.99) {
  tmp <- cov2cor(cov_mat)  
  if(any(tmp[lower.tri(tmp)] > limit)) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

