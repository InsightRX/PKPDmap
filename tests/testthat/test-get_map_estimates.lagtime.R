test_that("MAP fit works with oral model with lagtime", {
  ## Simulate some data
  reg <- PKPDsim::new_regimen(
    amt = 100,
    times = c(0, 24),
    type = "bolus"
  )
  data <- PKPDsim::sim(
    mod_1cmt_oral_lagtime,
    regimen = reg,
    parameters = list(CL = 5, V = 50, KA = 0.5, TLAG = 0.83),
    t_obs = c(0.25, 0.5, 0.75, 1.0, 2, 5, 8, 16),
    only_obs = TRUE
  )
  ## check that fit works for vector of lagtimes
  lagtimes <- c(0.25, 0.5, 0.75, 1.0, 2.0)
  est <- c()
  for(i in seq(lagtimes)) {
    fit <- get_map_estimates(
      model = mod_1cmt_oral_lagtime,
      parameters = list(CL = 6, V = 40, KA = 0.6, TLAG = lagtimes[i]), # slightly tweaked
      regimen = reg,
      omega = c(
        0.1, 
        0, 0.1,
        0, 0, 0.1
      ),
      fixed = c("KA"),
      error = list(prop = 0.1, add = 0.15),
      data = data,
      ## !important! due to lagtime discontinuity, cannot use gradient-based optimizers such as
      ## BGFS or L-BFGS-B, eventhough gradients are based on difference method. Optimization
      ## with these methods sometimes still work if lagtime initial estimate is higher than
      ## true value, but better to use a gradient-free method such as Nelder-Mead. 
      method = "Nelder-Mead",
      verbose = F
    )
    est <- c(est, fit$parameters$TLAG)
  }
  expect_equal(
    round(est, 3),
    c(0.594, 0.757, 0.856, 1.043, 1.307)
  )
})

test_that("np_hybrid fit works with explicit lagtime argument", {
  mod <- mod_1cmt_oral_lagtime
  reg <- PKPDsim::new_regimen(
    amt = 100,
    times = c(0, 24),
    type = "bolus"
  )
  lagtime <- 0.5
  data <- PKPDsim::sim(
    ode = mod,
    parameters = list(CL = 5, V = 50, KA = 0.5, TLAG = 0.85),
    regimen = reg,
    lagtime = lagtime,
    t_obs = c(1, 2, 4, 8, 12, 20, 25, 28, 32),
    only_obs = TRUE
  )
  fit <- get_map_estimates(
    model = mod,
    parameters = list(CL = 7, V = 60, KA = 0.7, TLAG = 0.85),
    regimen = reg,
    omega = c(0.05,
              0, 0.05,
              0, 0, 0.1),
    error = list(prop = 0.1, add = 0.1),
    data = data,
    lagtime = lagtime,
    fixed = "TLAG",
    type = "np_hybrid",
    verbose = FALSE
  )
  expect_true(!is.null(fit$np))
  expect_true(is.numeric(fit$parameters$CL))
  expect_true(is.numeric(fit$parameters$V))
})
