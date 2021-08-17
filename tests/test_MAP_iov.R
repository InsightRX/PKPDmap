library(testit)
library(PKPDmap)
library(PKPDsim)
Sys.setenv("R_TESTS" = "")

## Test IOV (in ODE block)
pars_iov   <- list(CL = 25, V = 55, CL_occ1 = 0, CL_occ2 = 0, CL_occ3 = 0, CL_occ4 = 0) 
model <- new_ode_model(code = "
                       CL_I = CL;
                       if(t<24) { CL_I = CL_I * exp(CL_occ1); }
                       if(t>=24) { if (t < 48) { CL_I = CL_I * exp(CL_occ2); } }
                       if(t>=48) { if (t < 72) { CL_I = CL_I * exp(CL_occ3); } }
                       if(t>=72) { CL_I = CL_I * exp(CL_occ4); }
                       dAdt[1] = -(CL_I/V) * A[1];
                       ", dose_code = NULL,
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
  as_eta = as_eta,
  verbose = T
)

## check that IOV etas moved away from initial estimate
testit::assert(round(fit$parameters$CL_occ1,2) == -0.15)
testit::assert(round(fit$parameters$CL_occ2,2) == 0.01)

## Test IOV (in PK block)
pars_iov   <- list(CL = 25, V = 55, CL_occ1 = 0, CL_occ2 = 0, CL_occ3 = 0, CL_occ4 = 0) 
model <- new_ode_model(code = "
                       dAdt[1] = -(CL_I/V) * A[1];
                       ", pk_code =   
                       "CL_I = CL;
                       if(times[i]<24) { CL_I = CL_I * exp(CL_occ1); }
                       if(times[i]>=24) { if (times[i] < 48) { CL_I = CL_I * exp(CL_occ2); } }
                       if(times[i]>=48) { if (times[i] < 72) { CL_I = CL_I * exp(CL_occ3); } }
                       if(times[i]>=72) { CL_I = CL_I * exp(CL_occ4); }
                       ",
                       declare_variables = c("CL_I"), parameters = pars_iov, obs = list(cmt = 1, scale = "V/1000"))

fit2 <- get_map_estimates(
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
