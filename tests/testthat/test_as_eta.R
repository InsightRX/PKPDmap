test_that("as_eta feature works", {
  model <- PKPDsim::new_ode_model(
    code = "
      CLi = CL * exp(et1);
      Vi = V * exp(et2);
      dAdt[1] = -(CLi/Vi) * A[1];", 
    obs = list(cmt = 1, scale = "V * exp(et2) /1000"),
    declare_variables = c("CLi", "Vi")
  )
  par <- list(
    CL = 15,
    V = 130,
    et1 = 0,
    et2 = 0
  )
  reg <- PKPDsim::new_regimen(
    amt = 100,
    n = 12,
    interval = 12,
    type = "infusion",
    t_inf = 1
  )
  obs <- PKPDsim::sim(
    ode = model,
    parameters = par,
    regimen = reg,
    only_obs = TRUE,
    t_obs = c(1, 11.5, 13, 23)
  )
  
  fit1 <- get_map_estimates(
    model = model,
    data = obs,
    parameters = list(
      CL = 10,
      V = 100,
      et1 = 0,
      et2 = 0
    ),
    fixed = c("et1", "et2"),
    regimen = reg,
    omega = c(0.1, 0.05, 0.1),
    error = list(prop = 0.1, add = 10),
    verbose = F
  )
  fit2 <- get_map_estimates(
    model = model,
    data = obs,
    parameters = list(
      CL = 10,
      V = 100,
      et1 = 0,
      et2 = 0
    ),
    fixed = c("CL", "V"),
    as_eta = c("et1", "et2"),
    regimen = reg,
    omega = c(0.1, 0.05, 0.1),
    error = list(prop = 0.1, add = 10),
    verbose = F
  )
  expect_equal(
    signif(fit1$parameters$CL, 4), signif(fit2$parameters$CL * exp(fit2$parameters$et1), 4)
  )
  expect_equal(
    signif(fit1$parameters$V, 4), signif(fit2$parameters$V * exp(fit2$parameters$et2), 4)
  )
})
