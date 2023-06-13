test_that("allow_obs_before_first_dose works", {
  mod <- PKPDsim::new_ode_model("pk_1cmt_iv")
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
