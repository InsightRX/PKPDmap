test_that("parse_weights returns default weights when weights is NULL", {
  t_obs <- c(0, 1, 2, 3)
  
  result <- parse_weights(NULL, t_obs)
  expect_equal(result, c(1, 1, 1, 1))
})

test_that("parse_weights returns default weights for single observation", {
  t_obs <- 0
  
  result <- parse_weights(NULL, t_obs)
  expect_equal(result, 1)
})

test_that("parse_weights returns default weights for empty observation vector", {
  t_obs <- numeric(0)
  
  result <- parse_weights(NULL, t_obs)
  expect_equal(result, numeric(0))
})

test_that("parse_weights returns provided weights when valid", {
  t_obs <- c(0, 1, 2, 3)
  weights <- c(1.5, 2.0, 0.5, 1.0)
  
  result <- parse_weights(weights, t_obs)
  expect_equal(result, weights)
})

test_that("parse_weights returns provided weights for single observation", {
  t_obs <- 0
  weights <- 2.5
  
  result <- parse_weights(weights, t_obs)
  expect_equal(result, 2.5)
})

test_that("parse_weights returns provided weights for empty observation vector", {
  t_obs <- numeric(0)
  weights <- numeric(0)
  
  result <- parse_weights(weights, t_obs)
  expect_equal(result, numeric(0))
})

test_that("parse_weights handles integer weights", {
  t_obs <- c(0, 1, 2)
  weights <- c(1L, 2L, 3L)
  
  result <- parse_weights(weights, t_obs)
  expect_equal(result, c(1, 2, 3))
})

test_that("parse_weights handles decimal weights", {
  t_obs <- c(0, 1, 2, 3)
  weights <- c(0.1, 0.5, 1.5, 2.0)
  
  result <- parse_weights(weights, t_obs)
  expect_equal(result, weights)
})

test_that("parse_weights handles zero weights", {
  t_obs <- c(0, 1, 2)
  weights <- c(0, 1, 0)
  
  result <- parse_weights(weights, t_obs)
  expect_equal(result, c(0, 1, 0))
})

test_that("parse_weights handles negative weights", {
  t_obs <- c(0, 1, 2)
  weights <- c(-1, 0, 1)
  
  result <- parse_weights(weights, t_obs)
  expect_equal(result, c(-1, 0, 1))
})

test_that("parse_weights handles large weights", {
  t_obs <- c(0, 1, 2)
  weights <- c(1000, 2000, 3000)
  
  result <- parse_weights(weights, t_obs)
  expect_equal(result, c(1000, 2000, 3000))
})

test_that("parse_weights handles very small weights", {
  t_obs <- c(0, 1, 2)
  weights <- c(1e-10, 1e-5, 1e-3)
  
  result <- parse_weights(weights, t_obs)
  expect_equal(result, c(1e-10, 1e-5, 1e-3))
})

test_that("parse_weights handles all equal weights", {
  t_obs <- c(0, 1, 2, 3, 4)
  weights <- rep(2.5, 5)
  
  result <- parse_weights(weights, t_obs)
  expect_equal(result, rep(2.5, 5))
})

test_that("parse_weights handles mixed positive and negative weights", {
  t_obs <- c(0, 1, 2, 3)
  weights <- c(1, -1, 0, 0.5)
  
  result <- parse_weights(weights, t_obs)
  expect_equal(result, c(1, -1, 0, 0.5))
})

test_that("parse_weights handles large observation vectors", {
  n <- 1000
  t_obs <- 1:n
  weights <- rep(1.5, n)
  
  result <- parse_weights(weights, t_obs)
  expect_equal(result, rep(1.5, n))
})

test_that("parse_weights handles large observation vectors with NULL weights", {
  n <- 1000
  t_obs <- 1:n
  
  result <- parse_weights(NULL, t_obs)
  expect_equal(result, rep(1, n))
})

test_that("parse_weights fails when weights vector is shorter than t_obs", {
  t_obs <- c(0, 1, 2, 3)
  weights <- c(1, 2)  # Only 2 weights for 4 observations
  
  expect_error(
    parse_weights(weights, t_obs),
    "Vector of weights of different size than observation vector!"
  )
})

test_that("parse_weights fails when weights vector is longer than t_obs", {
  t_obs <- c(0, 1)
  weights <- c(1, 2, 3, 4)  # 4 weights for 2 observations
  
  expect_error(
    parse_weights(weights, t_obs),
    "Vector of weights of different size than observation vector!"
  )
})

test_that("parse_weights fails when weights is empty but t_obs is not", {
  t_obs <- c(0, 1, 2)
  weights <- numeric(0)
  
  expect_error(
    parse_weights(weights, t_obs),
    "Vector of weights of different size than observation vector!"
  )
})

test_that("parse_weights fails when weights is not empty but t_obs is", {
  t_obs <- numeric(0)
  weights <- c(1, 2, 3)
  
  expect_error(
    parse_weights(weights, t_obs),
    "Vector of weights of different size than observation vector!"
  )
})

test_that("parse_weights handles NA values in weights", {
  t_obs <- c(0, 1, 2)
  weights <- c(1, NA, 3)
  
  result <- parse_weights(weights, t_obs)
  expect_equal(result, c(1, NA, 3))
})

test_that("parse_weights handles Inf values in weights", {
  t_obs <- c(0, 1, 2)
  weights <- c(1, Inf, 3)
  
  result <- parse_weights(weights, t_obs)
  expect_equal(result, c(1, Inf, 3))
})

test_that("parse_weights handles NaN values in weights", {
  t_obs <- c(0, 1, 2)
  weights <- c(1, NaN, 3)
  
  result <- parse_weights(weights, t_obs)
  expect_equal(result, c(1, NaN, 3))
})

test_that("parse_weights preserves weights data type", {
  t_obs <- c(0, 1, 2)
  
  # Test with numeric weights
  weights_numeric <- c(1.0, 2.0, 3.0)
  result_numeric <- parse_weights(weights_numeric, t_obs)
  expect_equal(result_numeric, weights_numeric)
  
  # Test with integer weights
  weights_integer <- c(1L, 2L, 3L)
  result_integer <- parse_weights(weights_integer, t_obs)
  expect_equal(result_integer, c(1, 2, 3))  # Note: integers become numeric
}) 

