mod <- PKPDsim::new_ode_model("pk_1cmt_iv")

test_that("allow_obs_before_first_dose works", {
  data <- data.frame(
    t = c(-0.45, 35.55),
    y = c(21.7, 17.8),
    evid = c(0, 0),
    loq = c(0, 0),
    obs_type = c(1, 1),
    unit = c("mg_l", "mg_l"),
    dv = c(21.7, 17.8),
    lastdose_time = c("2023-06-13 14:27:00", "2023-06-13 14:27:00"),
    id = c("64887d2aea22e60012333bb1", "64887d2aea22e60012333bb0"),
    ref_dose = c(NA, NA)
  )
  parameters <- structure(
    list(
      CL = 4.5,
      V = 58.4,
      V2 = 38.4,
      Q = 6.5,
      TH_CRCL = 0.8,
      TH_DIAL_CL = 0.7,
      TH_DIAL_V = 0.5,
      TDM_INIT = 21.7
    ),
    units = list(
      CL = "L/hr",
      V = "L/70kg",
      Q = "L/hr",
      V2 = "L",
      CLi = "L/hr",
      Vi = "L",
      Qi = "L/hr",
      V2i = "L"
    )
  )
  covariates <- list(
    WT = PKPDsim::new_covariate(value = 67.5, unit = "kg"),
    SEX = PKPDsim::new_covariate(value = 0),
    AGE = PKPDsim::new_covariate(value = 55.7, unit = "years"),
    CR = PKPDsim::new_covariate(value = c(0.67, 0.65), times = c(0, 23.3), unit = "mg_dl"),
    DIAL = PKPDsim::new_covariate(value = 0),
    CL_HEMO = PKPDsim::new_covariate(value = 0)
  )
  fixed <- c("Q", "TH_CRCL", "TH_DIAL_CL", "TH_DIAL_V", "TDM_INIT")
  omega <- c(0.1584, 0, 0.6659, 0, 0, 0.326)
  error <- list(prop = 0.227, add = 3.4)
  regimen <- structure(
    list(
      interval = 12,
      n = 4L,
      type = c("infusion", "infusion", "infusion", "infusion"),
      t_inf = c(1, 1, 1, 1),
      dose_times = c(0, 12, 24, 36),
      dose_amts = c(1000, 1000, 1000, 1000),
      first_dose_time = structure(
        1686666420,
        class = c("POSIXct", "POSIXt"),
        tzone = "UTC"
      )
    ),
    class = c("regimen", "list")
  )

  expect_error(
    get_map_estimates(
      model = mod,
      data = data,
      parameters = parameters,
      covariates = covariates,
      fixed = fixed,
      as_eta = NULL,
      omega = omega,
      error = error,
      obs_type_label = "obs_type",
      weights = c(0, 1),
      regimen = regimen,
      t_init = 0.45,
      allow_obs_before_dose = TRUE,
      A_init = c(0, 0, 0)
    ),
    NA
  )
})

test_that("MAP fit works with <LOQ data", {
  dat <- read.table(
    file = test_path("nm", "pktab1"),
    skip = 1,
    header = TRUE
  )
  colnames(dat)[1:3] <- c("id", "t", "y")
  data1 <- dat[dat$id == 1 & dat$EVID == 0, ]
  lloq <- 500
  data1$loq <- 0
  data1[data1$y < lloq, ]$loq <- -1
  data1[data1$loq < 0, ]$y  <- lloq
  reg <- PKPDsim::new_regimen(
    amt = 100000,
    times = c(0, 24),
    type = "bolus"
  )
  fit <- get_map_estimates(
    parameters = list(CL = 7.67, V = 97.7),
    model = mod,
    regimen = PKPDsim::new_regimen(
      amt = 100000,
      times = c(0, 24),
      type = "bolus"
    ),
    omega = c(0.0406, 0.0623, 0.117),
    error = list(prop = 0, add = sqrt(1.73E+04)),
    data = data1,
    censoring = "loq"
  )
  expect_equal(round(fit$parameters$CL, 3), 7.469)
  expect_equal(round(fit$parameters$V, 3), 92.382)
  expect_equal(
    round(fit$iwres, 3), 
    c(0.443, 1.656, -0.624, -0.91,
      -0.873,-0.579, 0.108, 0.754,
      0.095, 1.05, 0.122)
  )
  expect_equal(
    is.na(fit$wres),
    c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE,
      TRUE, FALSE, TRUE
    )
  )
})

test_that("MAP fit works with >ULOQ and <BLOQ data", {
  dat <- read.table(
    file = test_path("nm", "pktab1"),
    skip = 1,
    header = TRUE
  )
  colnames(dat)[1:3] <- c("id", "t", "y")
  data1 <- dat[dat$id == 1 & dat$EVID == 0, ]
  data1$loq <- 0
  uloq <- 1000
  lloq <- 500
  data1[data1$y < lloq, ]$loq <- -1
  data1[data1$loq < 0, ]$y  <- lloq
  data1[data1$y > uloq, ]$loq <- 1
  data1[data1$loq > 0, ]$y <- uloq
  fit <- get_map_estimates(
    parameters = list(CL = 7.67, V = 97.7),
    model = mod,
    regimen = PKPDsim::new_regimen(
      amt = 100000,
      times = c(0, 24),
      type = "bolus"
    ),
    omega = c(0.0406, 0.0623, 0.117),
    error = list(prop = 0, add = sqrt(1.73E+04)),
    data = data1,
    censoring = "loq"
  )
  expect_equal(round(fit$parameters$CL, 2), 7.54)
  expect_equal(round(fit$parameters$V, 1), 95.6)
  expect_equal(
    round(fit$iwres, 3),
    c(1.31, 1.204, -0.379, -0.698, -0.716,
      -0.465,  0.188, 0.735, 0.097, 1.068,
      0.126)
  )
  expect_equal(
    is.na(fit$wres), 
    c(TRUE, TRUE, FALSE, FALSE, FALSE, FALSE,
      FALSE, TRUE, TRUE, FALSE, TRUE
    )
  )
})

test_that("MAP fit works with LTBS res.error model", {
  dat <- read.table(
    file = test_path("nm", "pktab1"),
    skip = 1,
    header = TRUE
  )
  colnames(dat)[1:3] <- c("id", "t", "y")
  data1 <- dat[dat$id == 1 & dat$EVID == 0,]
  par   <- list(CL = 7.67, V = 97.7) 
  lloq <- 500
  data1$lloq <- 0
  data1[data1$y < lloq,]$lloq <- 1
  data1[data1$y < lloq,]$y <- lloq
  omega <- c(0.0406, 
             0.0623, 0.117)
  reg <- PKPDsim::new_regimen(amt = 100000, times=c(0, 24), type="bolus")
  fit1 <- get_map_estimates(
    parameters = par,
    model = mod,
    regimen = reg,
    omega = omega,
    error = list(prop = 0.15),
    data = data1,
    ltbs = FALSE
  )
  fit2 <- get_map_estimates(
    parameters = par,
    model = mod,
    regimen = reg,
    omega = omega,
    error = list(add = 0.15),
    data = data1,
    ltbs = TRUE
  )
  expect_equal(round(fit2$parameters$CL,3), 5.408)
  expect_equal(round(fit2$parameters$V,3), 98.954)  
})

test_that("MAP works for models with IOV", {
  pars_iov   <- list(
    CL = 25, V = 55, 
    CL_occ1 = 0, CL_occ2 = 0, CL_occ3 = 0, CL_occ4 = 0
  ) 
  model <- PKPDsim::new_ode_model(
    code = "
      CL_I = CL;
      if(t<24) { CL_I = CL_I * exp(CL_occ1); }
      if(t>=24) { if (t < 48) { CL_I = CL_I * exp(CL_occ2); } }
      if(t>=48) { if (t < 72) { CL_I = CL_I * exp(CL_occ3); } }
      if(t>=72) { CL_I = CL_I * exp(CL_occ4); }
      dAdt[1] = -(CL_I/V) * A[1];
      ", 
    dose_code = NULL,
    declare_variables = c("CL_I"), 
    parameters = pars_iov, 
    obs = list(cmt = 1, scale = "V/1000")
  )
  reg <- PKPDsim::new_regimen(
    amt = c(10.7, 10.7,10.7, 10.7, 10.7, 6.6,6.6,6.6),
    n = 8, 
    interval = 6, 
    type = 'infusion', 
    t_inf = 2
  )
  data <- data.frame(
    t = c(11.75, 14, 14.25, 14.5, 16, 18, 35.75, 38.1, 38.25, 40, 42),
    y = c(463, 1460, 1230, 1140, 714, 382, 284, 796, 544, 337, 222)/10,
    evid = 0
  )
  omega <- c(
    0.0531,
    0.0268, 0.0261, 
    0.0000, 0.0000, 0.1000, 
    0.0000, 0.0000, 0.0000, 0.1000
  )
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
  expect_equal(round(fit$parameters$CL_occ1,2), -0.15)
  expect_equal(round(fit$parameters$CL_occ2,2), 0.01)
  
  ## Test IOV (in PK block)
  pars_iov <- list(CL = 25, V = 55, CL_occ1 = 0, CL_occ2 = 0, CL_occ3 = 0, CL_occ4 = 0) 
  model <- PKPDsim::new_ode_model(
    code = "
      dAdt[1] = -(CL_I/V) * A[1];
    ", 
    pk_code = "
      CL_I = CL;
      if(times[i]<24) { CL_I = CL_I * exp(CL_occ1); }
      if(times[i]>=24) { if (times[i] < 48) { CL_I = CL_I * exp(CL_occ2); } }
      if(times[i]>=48) { if (times[i] < 72) { CL_I = CL_I * exp(CL_occ3); } }
       if(times[i]>=72) { CL_I = CL_I * exp(CL_occ4); }
    ",
    declare_variables = c("CL_I"), 
    parameters = pars_iov, 
    obs = list(cmt = 1, scale = "V/1000")
  )
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
  expect_equal(round(fit$parameters$CL_occ1,2), -0.15)
  expect_equal(round(fit$parameters$CL_occ2,2), 0.01)
})

test_that("Test if datasets with duplicate measurements work.", {
  ## Rationale:
  ## - users might want to use this to put more weight on samples (preferrable is to use weighting functionality, obviously)
  ## - data (e.g. from EMR) might have duplicates 
  ## We should show a warning that there are duplicates, but the estimation should not fail.
  dat <- read.table(
    file = test_path("nm", "pktab1"),
    skip = 1,
    header = TRUE
  )
  colnames(dat)[1:3] <- c("id", "t", "y")
  tmp   <- dat[dat$id == 1 & dat$EVID == 0,]
  par   <- list(CL = 7.67, V = 97.7) 
  fits <- c()
  omega <- c(0.0406, 
             0.0623, 0.117)
  reg <- PKPDsim::new_regimen(amt = 100000, times=c(0, 24), type="bolus")
  
  ## add duplicates at t=4 and t=24
  ## make it difficult, becasue we'll put t=24 exactly at same timepoint as new dose
  ##
  tmp[tmp$t == 23.9,]$t <- 24 
  tmp <- rbind(tmp, tmp[tmp$t %in% c(4, 24),])
  tmp <- tmp[order(tmp$t),]
  
  expect_message(
    fit <- get_map_estimates(
      parameters = par,
      model = mod,
      regimen = reg,
      omega = omega,
      weights = rep(1, length(tmp$t)),
      error = list(prop = 0, add = sqrt(1.73E+04)),
      data = tmp
    )
  )
  expect_true("map_estimates" %in% class(fit))  
})

test_that("When vcov matrix is not postitive definite, don't crash but return original omega matrix and throw warning", {
  ## Sometimes the Hessian can't be inverted and the vcov matrix is then not positive definite.
  ## In that case get_map_estimates() should just return the original omega matrix,
  ## and throw a warning.
  
  ## The below code creates a TDM scenario that will result in such a situation,
  ## this is with the vancomycin Carreno model and a single TDM.
  par <- list(
    V = 25.76,
    SCLSlope = 0.036,
    K12 = 2.29,
    K21 = 1.44,
    SCLInter = 0.18,
    TDM_INIT = 0
  )
  omega <- c(0.205590,
             0.000000, 0.308640,
             0.000000, 0.000000, 1.116760,
             0.000000, 0.000000, 0.000000, 1.443300,
             0.000000, 0.000000, 0.000000, 0.000000, 0.057068)
  
  ## Set up model and regimen
  # egfr <- clinPK::calc_egfr(age = 65, weight = 156, height = 180, scr = 1.49, sex = "male", method = "cockcroft_gault_adjusted")$value
  egfr <- 75.080589758495
  covs <- list(CRCL = PKPDsim::new_covariate(value = egfr))
  model <- PKPDsim::new_ode_model(
    code = "
    CLi = SCLInter + SCLSlope * (CRCL*16.667) \
    Vi = V \
    Qi = K12 * Vi \
    V2i = Qi / K21 \
    dAdt[0] = -(CLi/V)*A[0] - K12*A[0] + K21*A[1] \
    dAdt[1] = K12*A[0] - K21*A[1] \
    dAdt[2] = A[0]/V
  ",
    pk_code = "",
    parameters = par,
    omega_matrix = omega,
    obs = list(cmt = 1, scale="V"),
    dose = list(cmt = 1, bioav = 1),
    declare_variables = c("CLi", "Vi", "Qi", "V2i"),
    covariates = covs
  )
  reg <- PKPDsim::new_regimen(
    amt = c(1000, 1000),
    times = c(0, 19.75),
    t_inf = 1.5,
    type = "infusion")
  t_tdm <- 19.75 + 13.75
  data <- data.frame(t = t_tdm, y = 8.8)
  
  expect_warning(
    fit <- get_map_estimates(
      model = model,
      parameters = par,
      covariates = covs,
      regimen = reg,
      data = data,
      omega = omega,
      error = list(prop = 0.1, add = 0.1),
      fixed = "TDM_INIT"
    )
  )
  expect_equal(fit$vcov, omega)
})

test_that("Can estimate mixture model", {
  model <- PKPDsim::new_ode_model(
    code = "
      dAdt[0] = -(CL/V)*A[0];
    ",
    pk_code = " ",
    obs = list(cmt = 1, scale = "V"),
    mixture = list(CL = list(values = c(5, 9), probability = 0.3))
  )
  
  dat   <- read.table(file=test_path("nm", "pktab1"), skip=1, header=TRUE)
  colnames(dat)[1:3] <- c("id", "t", "y")
  data1 <- dat[dat$id == 1 & dat$EVID == 0,]
  par   <- list(CL = 7.67, V = 97.7) 
  omega <- c(0.0406, 
             0.0623, 0.117)
  fit <- get_map_estimates(parameters = par,
                           model = model,
                           regimen = PKPDsim::new_regimen(amt = 100000, times=c(0, 24), type="bolus"),
                           omega = omega,
                           error = list(prop = 0, add = sqrt(1.73E+04)),
                           data = data1)
  
  expect_equal(round(fit$parameters$CL,3), 5.368)
  expect_equal(round(fit$parameters$V,3), 99.762)
  expect_true(!is.null(fit$mixture$probabilities))
  expect_equal(round(fit$mixture$probabilities, 2), c(0.66, 0.34))
  expect_equal(fit$mixture$mixture_group, 1)
  expect_equal(fit$mixture$selected, 5)  
})

test_that("MAP fit works with TDM before first dose", {
  ## TDM before first dose:
  ## at t=-8, conc=10000
  ## Use this as true value in the simulations
  
  dat   <- read.table(file=test_path("nm", "pktab1"), skip=1, header=TRUE)
  colnames(dat)[1:3] <- c("id", "t", "y")
  dat   <- dat[dat$id == 1 & dat$EVID == 0,]
  par   <- list(CL = 7.67, V = 97.7, TDM_INIT = 500) 
  reg   <- PKPDsim::new_regimen(amt = 100000, times=c(0, 24), type="bolus")
  model <- PKPDsim::new_ode_model(
    code = "
      dAdt[1] = -(CL/V)*A[1];
    ", 
    state_init = "
      A[1] = TDM_INIT * V
    ",
    parameters = par, 
    obs = list(cmt = 1, scale = "V"), 
    cpp_show_code = F
  )
  fits <- c()
  omega <- c(0.0406, 
             0.0623, 0.117)
  fit4 <- get_map_estimates(parameters = par,
                            model = model,
                            fixed = c("TDM_INIT"),
                            regimen = reg,
                            omega = omega,
                            weights = rep(1, length(dat$t)),
                            error = list(prop = 0, add = sqrt(1.73E+04)),
                            data = dat,
                            t_init = 4,
                            verbose = F
  )
  
  fit2 <- get_map_estimates(parameters = par,
                            model = model,
                            fixed = c("TDM_INIT"),
                            regimen = reg,
                            omega = omega,
                            weights = rep(1, length(dat$t)),
                            error = list(prop = 0, add = sqrt(1.73E+04)),
                            data = dat,
                            t_init = 2,
                            verbose = F
  )
  fit0 <- get_map_estimates(parameters = par,
                            model = model,
                            fixed = c("TDM_INIT"),
                            regimen = reg,
                            omega = omega,
                            weights = rep(1, length(dat$t)),
                            error = list(prop = 0, add = sqrt(1.73E+04)),
                            data = dat,
                            verbose = F
  )
  expect_true("map_estimates" %in% class(fit4))
  expect_true("map_estimates" %in% class(fit2))
  expect_true("map_estimates" %in% class(fit0))
  expect_true(fit4$parameters$CL != fit2$parameters$CL)
  expect_true(fit2$parameters$CL != fit0$parameters$CL)
})


test_that("MAP works for multiple observation types", {
  ## define parameters
  pk1 <- PKPDsim::new_ode_model(code = "dAdt[1] = -(CL/V)*A[1]", obs = list(scale="V/1000", cmt=1))
  regimen  <- PKPDsim::new_regimen(amt = 100, interval = 12, n = 5, type="infusion", t_inf = 1)
  parameters   <- list("CL" = 15, "V" = 150)
  omega <- PKPDsim::cv_to_omega(list("CL" = 0.2, "V" = 0.2), parameters[1:2])
  
  ruv_single <- list(prop = 0.1, add = 1)
  ruv_multi <- list(prop = c(0.1, 1), add = c(0.1, 20))
  ruv_multi2 <- list(prop = c(0.1, 1, 2), add = c(0.1, 1, 20))
  
  ## simulate single individual in population
  # some observations with much higher residual error, should affect fit that much
  set.seed(83475)
  data_multi <- PKPDsim::sim_ode(
    ode = pk1, 
    parameters = list(CL = 20, V = 200), 
    regimen = regimen,
    int_step_size = 0.1,
    only_obs = TRUE,
    obs_type =  c(1,2,1,2), 
    t_obs = c(2, 4, 6, 8),
    res_var = ruv_multi
  )
  
  data_noruv <- PKPDsim::sim_ode(
    ode = pk1, 
    parameters = list(CL = 20, V = 200), 
    regimen = regimen,
    int_step_size = 0.1,
    only_obs = TRUE,
    t_obs = c(2, 4, 6, 8)
  )
  
  ## Fit
  # first fit with single low residual errror
  fit1 <- get_map_estimates(
    model = pk1, 
    data = data_multi, 
    parameters = parameters,
    omega = omega,
    regimen = regimen,
    obs_type_label = NULL,
    error = ruv_single
  )
  
  # then fit with multiple residual error. Fit should be affected less by points with large error
  fit2 <- get_map_estimates(
    model = pk1, 
    data = data_multi, 
    parameters = parameters,
    omega = omega,
    regimen = regimen,
    obs_type_label = "obs_type",
    error = ruv_multi
  )
  fit3 <- get_map_estimates(
    model = pk1, 
    data = data_noruv, 
    parameters = parameters,
    omega = omega,
    regimen = regimen,
    error = ruv_single
  )
  expect_true((abs(fit3$parameters$CL - fit1$parameters$CL) / abs(fit3$parameters$CL - fit2$parameters$CL)) > 5)
  
  ## Test ordering of data
  data_multi2 <- PKPDsim::sim_ode(
    ode = pk1, 
    parameters = list(CL = 20, V = 200), 
    regimen = regimen,
    int_step_size = 0.1,
    only_obs = TRUE,
    obs_type =  c(1, 2, 2, 1, 1, 2), 
    t_obs = c(2, 4, 6, 6, 8, 8),
    res_var = ruv_multi
  )
  data_multi3 <- PKPDsim::sim_ode(
    ode = pk1, 
    parameters = list(CL = 20, V = 200), 
    regimen = regimen,
    int_step_size = 0.1,
    only_obs = TRUE,
    obs_type =  c(1, 2, 3, 2, 1, 1, 2), 
    t_obs = c(2, 2, 2, 4, 4, 8, 8),
    res_var = ruv_multi2
  )
  
  # PKPDsim is expected to re-order. Testing here to make sure this is still the case
  expect_equal(data_multi2$obs_type, c(1, 2, 1, 2, 1, 2))
  expect_equal(data_multi3$obs_type, c(1, 2, 3, 1, 2, 1, 2))
  fit4 <- get_map_estimates(
    model = pk1, 
    data = data_multi3[rev(seq_len(nrow(data_multi3))), ], 
    parameters = parameters,
    omega = omega,
    regimen = regimen,
    obs_type_label = "obs_type",
    error = ruv_multi2
  )
  expect_equal(fit4$obs_type, data_multi3$obs_type)
})


