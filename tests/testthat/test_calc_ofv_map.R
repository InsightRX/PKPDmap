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
    weight_prior = 1, 
    include_omega = TRUE, 
    include_error = TRUE)
  expect_equal(round(res, 5), 
               c(-0.8171, 43.87438, 3.594, 9.28779))
})

## check that weights are applied
test_that("checkt that tdm weights are applied", {
  res2 <- calc_ofv_map(
    eta, 
    omega,
    omega_inv,
    omega_eigen,
    dv, 
    ipred, 
    res_sd,
    weights = c(0.2, 0.5, 1),
    weight_prior = 1, 
    include_omega = TRUE, 
    include_error = TRUE)
  expect_equal(round(res2, 5), 
               c(-0.8171, 8.77488, 1.797, 9.28779))
  
})

test_that("check that weight prior is applied correctly", {
  res3 <- calc_ofv_map(
    eta, 
    omega,
    omega_inv,
    omega_eigen,
    dv, 
    ipred, 
    res_sd,
    weights = 0.3,
    weight_prior = 1, 
    include_omega = TRUE, 
    include_error = TRUE)
  expect_equal(round(res3, 5), 
               c(-0.8171, 13.16231, 1.0782, 2.78634))
})
