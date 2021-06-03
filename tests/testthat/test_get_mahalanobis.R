test_that("Mahalanobis is calculated - w_ipred length 1", {
  expect_equal(get_mahalanobis(1.9, 1.459, 0.478, FALSE), 0.85117995)
  expect_equal(get_mahalanobis(1.9, 1.459, 0.478, TRUE),  0.30527386)
})

test_that("Mahalanobis is calculated - w_ipred length >1", {
  y <- c(3500, 2800, 2200, 1400, 500)
  ipred <- c(3388, 2844, 2388, 1186, 589)
  w_ipred <- c(359, 302, 254, 127, 66)
  expect_equal(get_mahalanobis(y, ipred, w_ipred, FALSE), 5.3241594)
  expect_equal(get_mahalanobis(y, ipred, w_ipred, TRUE),  7.982002e-06)
})

test_that("Computationally singular cov matrix returns NULL", {
  y <- c(0, 0, 0, 50, 27.5)
  ipred <- c(2.16e-10, 1.70e-09, 7.96e-08, -993, -10611099252899028)
  w_ipred <- c(2.81e-11, 2.21e-10, 1.03e-08, 1.29e+02, 1.38e+15)
  expect_null(get_mahalanobis(y, ipred, w_ipred, FALSE))
})
