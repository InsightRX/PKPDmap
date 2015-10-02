library(testit)
library(PKPDmap)
library(PKPDsim)
Sys.setenv("R_TESTS" = "")

## Basic precision and accuracy of MAP estimation (compared to NONMEM)
dat   <- read.table(file=paste0(system.file(package="PKPDmap"), "/pktab1"), skip=1, header=TRUE)
ebe   <- read.csv(file=paste0(system.file(package="PKPDmap"), "/patab1.csv"))
colnames(dat)[1:3] <- c("id", "t", "y")
dat   <- dat[dat$id <= 20,] # not necessary to do all 100 id's, if bug will be clear from 20 ids too.
ebe   <- ebe[ebe$id <= 20,]  
par   <- list(CL = 7.67, V = 97.7) 

model <- new_ode_model("pk_1cmt_iv")
fits <- c()
for(i in seq(unique(dat$id))) {
  tmp <- get_map_estimates(parameters = par,
                           model = model,
                           regimen = new_regimen(amt = 100000, times=c(0, 24), type="bolus"),
                           omega = c(0.0406, 0.0623, 0.117),
                           error = list(prop = 0, add = sqrt(1.73E+04)),
                           data = dat[dat$id == i,])
  fits <- rbind(fits, cbind(tmp$parameters$CL, tmp$parameters$V))
}

## allow max deviation of 0.5% with NONMEM
assert("PK 1cmt iv MAP CL estimate same as NONMEM (delta < 0.5%)", 
  max(fits[,1] - ebe$CL) / mean(ebe$CL) < 0.005)
assert("PK 1cmt iv MAP V estimate same as NONMEM (delta < 0.5%)", 
  max(fits[,2] - ebe$V) / mean(ebe$V) < 0.005)


