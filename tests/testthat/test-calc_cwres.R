test_that("calc_cwres FOCE vcov uses per-observation weights", {
  # Minimal 1-compartment IV model with two observations.
  model <- PKPDsim::new_ode_model(
    code = "dAdt[0] = -(CL/V) * A[0]",
    obs = list(cmt = 1, scale = "V"),
    dose = list(cmt = 1, bioav = 1),
    parameters = list(CL = 5, V = 50)
  )
  regimen <- PKPDsim::new_regimen(
    amt = 1000, times = 0, t_inf = 1, type = "infusion"
  )

  parameters_population <- list(CL = 5, V = 50)
  nonfixed <- c("CL", "V")
  omega_full <- matrix(c(0.1, 0, 0, 0.1), nrow = 2)
  eta_hat <- c(0, 0)
  t_obs <- c(2, 8)

  # Individual predictions at eta_hat (population params, since eta = 0).
  ipred_raw <- {
    suppressMessages(
      PKPDsim::sim_ode(
        ode = model, parameters = parameters_population, n_ind = 1,
        regimen = regimen, t_obs = t_obs, only_obs = TRUE, checks = FALSE
      )$y
    )
  }
  y <- ipred_raw  # observations equal to predictions (residuals irrelevant for vcov)

  args <- list(
    eta_hat = eta_hat,
    ipred_raw = ipred_raw,
    y = y,
    omega_full = omega_full,
    error = list(prop = 0.1, add = 0.1),
    obs_type = c(1, 1),
    transf = function(x) x,
    model = model,
    parameters_population = parameters_population,
    nonfixed = nonfixed,
    as_eta = c(),
    covariates = NULL,
    regimen = regimen,
    lagtime = NULL,
    t_obs = t_obs,
    obs_type_sim = c(1, 1)
  )

  res_unit <- do.call(calc_cwres, c(args, list(weights = rep(1, 2))))
  res_weighted <- do.call(calc_cwres, c(args, list(weights = rep(2, 2))))

  # Higher observation weights increase the likelihood term in the FOCE
  # Hessian H = Omega^-1 + F' diag(w/sigma^2) F, so vcov = H^-1 shrinks.
  # Without weights threaded into the Hessian the two vcovs would be identical.
  expect_false(is.null(res_unit$vcov))
  expect_false(is.null(res_weighted$vcov))
  expect_false(isTRUE(all.equal(res_unit$vcov, res_weighted$vcov)))
  expect_true(det(res_weighted$vcov) < det(res_unit$vcov))

  # CWRES (a model diagnostic) always uses unit weights, so it is unchanged.
  expect_equal(res_unit$cwres, res_weighted$cwres)
})
