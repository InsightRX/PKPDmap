eta <- matrix(c(0.1, -0.1), nrow=1)
omega <- matrix(c(0.1, 0.05, 
                  0.05, 0.1), ncol=2)
omega_inv <- solve(omega)
omega_eigen <- sum(log(eigen(omega)$values))
dv <- c(15, 20, 25)
ipred <- c(17, 19, 23)
res_sd <- (dv-ipred) * 0.1 + 0.5

test_that("Basic OFV calc works", {
  res <- calc_ofv_map(
    eta, 
    omega,
    omega_inv,
    omega_eigen,
    dv, 
    ipred, 
    res_sd,
    weights = 1,
    include_omega = TRUE, 
    include_error = TRUE)
  expect_equal(round(res, 5), 
               c(-0.8171, 43.87438, 3.594, 9.28779))
})

test_that("check that tdm weights are applied", {
  ## This is weights of observations, not prior weight
  ## Prior weight is not an argument to calc_ofv_map(), this is handled in get_map_estimates()
  res2a <- calc_ofv_map(
    eta, 
    omega,
    omega_inv,
    omega_eigen,
    dv, 
    ipred, 
    res_sd,
    weights = 1,
    include_omega = TRUE, 
    include_error = TRUE
  )
  res2b <- calc_ofv_map(
    eta, 
    omega,
    omega_inv,
    omega_eigen,
    dv, 
    ipred, 
    res_sd,
    weights = c(0.2, 0.5, .8),
    include_omega = TRUE, 
    include_error = TRUE
  )
  expect_true(
    all(round(res2a[-1], 5) !=  round(res2b[-1], 5))
  )
})

test_that("check that include_omega works", {
  res4a <- calc_ofv_map(
    eta, 
    omega,
    omega_inv,
    omega_eigen,
    dv, 
    ipred, 
    res_sd,
    weights = 1,
    include_omega = FALSE, 
    include_error = TRUE)
  res4b <- calc_ofv_map(
    eta, 
    omega * 2,
    omega_inv,
    omega_eigen,
    dv, 
    ipred, 
    res_sd,
    weights = 1,
    include_omega = FALSE, 
    include_error = TRUE)
  expect_equal(res4a, res4b)
})

test_that("check that include_error gives correct likelihoods", {
  res5 <- calc_ofv_map(
    eta, 
    omega,
    omega_inv,
    omega_eigen,
    dv, 
    ipred, 
    res_sd,
    weights = 1,
    include_omega = TRUE, 
    include_error = TRUE)
  expect_equal(round(res5, 5), c(-0.8171, 43.87438, 3.594, 9.28779))
})
