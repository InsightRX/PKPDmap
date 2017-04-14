library(testit)
library(PKPDmap)
library(PKPDsim)
Sys.setenv("R_TESTS" = "")

## Basic precision and accuracy of MAP estimation (compared to NONMEM)
dat   <- read.table(file=paste0(system.file(package="PKPDmap"), "/pktab1"), skip=1, header=TRUE)
ebe   <- read.csv(file=paste0(system.file(package="PKPDmap"), "/patab1.csv"))
sdtab <- read.table(file=paste0(system.file(package="PKPDmap"), "/sdtab1"), skip=1, header=TRUE)
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
                           omega = c(0.0406, 
                                     0.0623, 0.117),
                           weights = rep(1, length(dat[dat$id == i & dat$EVID == 0,1])),
                           error = list(prop = 0, add = sqrt(1.73E+04)),
                           data = dat[dat$id == i & dat$EVID == 0,]
                           #residuals = TRUE, verbose = TRUE
  )
  fits <- rbind(fits, cbind(tmp$parameters$CL, tmp$parameters$V))
}

## allow max deviation of 0.5% with NONMEM
assert("PK 1cmt iv MAP CL estimate same as NONMEM (delta < 0.5%)", 
       max(fits[,1] - ebe$CL) / mean(ebe$CL) < 0.005)
assert("PK 1cmt iv MAP V estimate same as NONMEM (delta < 0.5%)", 
       max(fits[,2] - ebe$V) / mean(ebe$V) < 0.005)

## Test 2cmt model (just check if it runs, no asserts)
model2 <- new_ode_model("pk_2cmt_iv")
fits <- c()
par2 <- list(CL = 7.67, V = 97.7, Q = 3, V2 = 50)
for(i in seq(unique(dat$id))) {
  tmp <- get_map_estimates(parameters = par2,
                           model = model2,
                           regimen = new_regimen(amt = 100000, times=c(0, 24), type="bolus"),
                           omega = c(0.0406,
                                     0.0623, 0.117,
                                     0.001, 0.01, 0.1,
                                     0.001, 0.001, 0.01, 0.1),
                           error = list(prop = 0, add = sqrt(1.73E+04)),
                           int_step_size = 0.1,
                           data = dat[dat$id == i,])
  fits <- rbind(fits, cbind(tmp$parameters$CL, tmp$parameters$V))
}

## Test 2cmt model but fix V2
par2 <- list(CL = 7.67, V = 97.7, Q = 3, V2 = 50)
i <- 1
tmp1 <- get_map_estimates(parameters = par2,
                          model = model2,
                          regimen = new_regimen(amt = 100000, times=c(0, 24), type="bolus"),
                          omega = c(0.0406,
                                    0.0623, 0.117,
                                    0.001, 0.01, 0.1,
                                    0.001, 0.01, 0.01, 0.1),
                          fixed = c("V2"),
                          error = list(prop = 0, add = sqrt(1.73E+04)),
                          int_step_size = 0.1,
                          data = dat[dat$id == i,])
tmp2 <- get_map_estimates(parameters = par2,
                          model = model2,
                          regimen = new_regimen(amt = 100000, times=c(0, 24), type="bolus"),
                          omega = c(0.0406,
                                    0.0623, 0.117,
                                    0.001, 0.01, 0.1),
                          fixed = c("V2"),
                          error = list(prop = 0, add = sqrt(1.73E+04)),
                          int_step_size = 0.1,
                          data = dat[dat$id == i,])
tmp3 <- get_map_estimates(parameters = par2,
                          model = model2,
                          regimen = new_regimen(amt = 100000, times=c(0, 24), type="bolus"),
                          omega = c(0.0406,
                                    0.0623, 0.117),
                          fixed = c("Q", "V2"),
                          error = list(prop = 0, add = sqrt(1.73E+04)),
                          int_step_size = 0.1,
                          data = dat[dat$id == i,])
assert("check fixing parameters #1",
       tmp1$parameters$V2 == 50)
assert("check fixing parameters #2",
       tmp2$parameters$V2 == 50 && tmp2$parameters$Q != 3)
assert("check fixing parameters #2",
       tmp3$parameters$V2 == 50 && tmp3$parameters$Q == 3)

assert("check fixing parameters #2",
       tmp3$parameters$V2 == 50 && tmp3$parameters$Q == 3)

tmp4 <- get_map_estimates(parameters = par2,
                          model = model2,
                          regimen = new_regimen(amt = 100000, times=c(0, 24), type="bolus"),
                          omega = c(0.0406,
                                    0.0623, 0.117),
                          fixed = c("Q", "V2"),
                          weights = weight_by_time(
                            dat[dat$id == i & dat$EVID == 0,]$t,
                            t_end_gradient = 12),
                          error = list(prop = 0, add = sqrt(1.73E+04)),
                          int_step_size = 0.1,
                          data = dat[dat$id == i,])

assert("weighting gives different results",
       tmp3$parameters$CL != tmp4$parameters
)
tmp5 <- get_map_estimates(parameters = par2,
                          model = model2,
                          regimen = new_regimen(amt = 100000, times=c(0, 24), type="bolus"),
                          omega = c(0.0406,
                                    0.0623, 0.117),
                          fixed = c("Q", "V2"),
                          weights = weight_by_time(
                            dat[dat$id == i & dat$EVID == 0,]$t,
                            t_end_gradient = 0),
                          error = list(prop = 0, add = sqrt(1.73E+04)),
                          int_step_size = 0.1,
                          data = dat[dat$id == i,])
assert("weighting with all weight=1 gives same results",
       tmp5$parameters$CL == tmp3$parameters$CL
)

grad <- weight_by_time(time = c(0:24), 
                       t_start_gradient = 5,
                       t_end_gradient = 23)
assert("gradient correct",
       sum(grad[1:6]) == 0 && sum(tail(grad,2)) == 2
)
assert("gradient correct",
       sum(weight_by_time(0:9)) == 5
)

## check individual residuals
id <- 1
tmp <- get_map_estimates(parameters = par,
                         model = model,
                         regimen = new_regimen(amt = 100000, times=c(0, 24), type="bolus"),
                         omega = c(0.0406, 
                                   0.0623, 0.117),
                         weights = rep(1, length(dat[dat$id == id & dat$EVID == 0,1])),
                         error = list(prop = 0, add = sqrt(1.73e04)),
                         data = dat[dat$id == id,], 
                         residuals = TRUE)
ires1 <- sdtab[sdtab$ID == 1 & sdtab$EVID == 0,]$IRES
iwres1 <- sdtab[sdtab$ID == 1 & sdtab$EVID == 0,]$IWRES
assert("residuals <1% different from reference", 
       all((tmp$residuals - ires1) / ires1 < 0.01))
assert("weighted residuals 0.01 different from reference", 
       all(abs(tmp$weighted_residuals - iwres1) < 0.01))


## Test IOV
pars_iov   <- list(CL = 25, V = 55, CL_occ1 = 0, CL_occ2 = 0, CL_occ3 = 0, CL_occ4 = 0) 
model <- new_ode_model(code = "
  dAdt[1] = -(CL_I/V) * A[1];
", dose_code = "
  CL_I = CL;
  if(times[i]<24) { CL_I = CL_I * exp(CL_occ1); }
  if(times[i]>=24) { if (times[i] < 48) { CL_I = CL_I * exp(CL_occ2); } }
  if(times[i]>=48) { if (times[i] < 72) { CL_I = CL_I * exp(CL_occ3); } }
  if(times[i]>=72) { CL_I = CL_I * exp(CL_occ4); }
",
declare_variables = c("CL_I"), parameters = pars_iov, obs = list(cmt = 1, scale = "V/1000"))

reg <- new_regimen(amt = c(10.7, 10.7,10.7, 10.7, 10.7, 6.6,6.6,6.6), n = 8, interval = 6, type = 'infusion', t_inf = 2)
data <- data.frame(cbind(t = c(11.75, 14, 14.25, 14.5, 16, 18, 
                               35.75, 38.1, 38.25, 40, 42),
                         y = c(463, 1460, 1230, 1140, 714, 382, 
                               284, 796, 544, 337, 222)/10,
                         evid = 0))
omega <- c(0.0531, 0.0268, 0.0261, 0.0000, 0.0000, 0.1000, 0.0000, 0.0000, 0.0000, 0.1000)
fixed <- c("CL_occ3", "CL_occ4")
ruv <- list(prop = 0.115, add = 21)
# sim(ode = model, parameters = pars_iov, regimen = reg)
as_eta <- c("CL_occ1", "CL_occ2", "CL_occ3", "CL_occ4")

fit <- get_map_estimates(
  model = model, 
  data = data,
  method = 'BFGS',
  omega = omega,
  parameters = pars_iov, 
  regimen = reg, 
  fixed = fixed,
  error = ruv,
  as_eta = as_eta
)

## check that IOV etas moved away from initial estimate
testit::assert(round(fit$parameters$CL_occ1,2) == -0.15)
testit::assert(round(fit$parameters$CL_occ2,2) == 0.01)

# library(PKPDplot)
# ipred <- sim(ode = model, parameters = fit$parameters, regimen = reg)
# ipred %>% plot() + geom_point(data = data, aes(x=t, y=y))


## Check that weighting works
model <- new_ode_model(code = "
  dAdt[1] = -(CL/V) * A[1];
             ", obs = list(cmt = 1, scale = "V/1000"))
par <- list(CL = 10, V = 100)
reg <- new_regimen(amt = 100, n = 12, interval = 12, type="infusion", t_inf = 1)
obs <- PKPDsim::sim(ode = model, parameters = par, regimen = reg, 
                    only_obs = TRUE, t_obs = c(11.5, 143.5))

obs$evid <- 0
obs$y[1] <- 600
fit1 <- get_map_estimates(model = model, data = obs, 
                          parameters = list(CL = 11, V = 90),
                          regimen = reg, omega = c(0.1, 0.05, 0.1),
                          error = list(prop = 0.1, add = 10), weights = c(1, 1), 
                          residuals = T)
assert("CL estimate no weighting", round(fit1$parameters$CL,1) == 6.4)

fit2 <- get_map_estimates(model = model, data = obs, 
                          parameters = list(CL = 11, V = 90),
                          regimen = reg, omega = c(0.1, 0.05, 0.1),
                          error = list(prop = 0.1, add = 10), weights = c(0.25, 1), 
                          residuals = T)
assert("CL estimate changed appropriately", round(fit2$parameters$CL,1) == 7.7)
assert("Scaling of residuals based on weighting", diff(abs(fit2$iwres)) < 0.1) # should be small difference. Scaling of residuals is not 100% correct, but seems close enough 
