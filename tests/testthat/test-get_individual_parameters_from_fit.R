test_that("get_individual_parameters_from_fit handles basic case", {
  # Create test data
  fit <- list(coef = c(0.1, 0.2))  # Two coefficients
  parameters <- list(CL = 10, V = 20)
  nonfixed <- c("CL", "V")
  as_eta <- NULL
  
  result <- get_individual_parameters_from_fit(fit, parameters, nonfixed, as_eta)
  
  # Check that parameters are correctly transformed
  expect_equal(result$CL, 10 * exp(0.1))
  expect_equal(result$V, 20 * exp(0.2))
})

test_that("get_individual_parameters_from_fit handles as_eta parameters", {
  # Create test data with some parameters as direct etas
  fit <- list(coef = c(0.1, 0.2))
  parameters <- list(CL = 10, V = 20)
  nonfixed <- c("CL", "V")
  as_eta <- c("CL")  # CL is direct eta
  
  result <- get_individual_parameters_from_fit(fit, parameters, nonfixed, as_eta)
  
  # Check that CL is direct eta and V is transformed
  expect_equal(result$CL, 0.1)
  expect_equal(result$V, 20 * exp(0.2))
})

test_that("get_individual_parameters_from_fit handles all as_eta parameters", {
  # Create test data with all parameters as direct etas
  fit <- list(coef = c(0.1, 0.2))
  parameters <- list(CL = 10, V = 20)
  nonfixed <- c("CL", "V")
  as_eta <- c("CL", "V")
  
  result <- get_individual_parameters_from_fit(fit, parameters, nonfixed, as_eta)
  
  # Check that all parameters are direct etas
  expect_equal(result$CL, 0.1)
  expect_equal(result$V, 0.2)
})

test_that("get_individual_parameters_from_fit handles no nonfixed parameters", {
  # Create test data with no nonfixed parameters
  fit <- list(coef = numeric(0))
  parameters <- list(CL = 10, V = 20)
  nonfixed <- character(0)
  as_eta <- NULL
  
  result <- get_individual_parameters_from_fit(fit, parameters, nonfixed, as_eta)
  
  # Check that parameters remain unchanged
  expect_equal(result$CL, 10)
  expect_equal(result$V, 20)
})

test_that("get_individual_parameters_from_fit handles negative coefficients", {
  # Create test data with negative coefficients
  fit <- list(coef = c(-0.1, -0.2))
  parameters <- list(CL = 10, V = 20)
  nonfixed <- c("CL", "V")
  as_eta <- NULL
  
  result <- get_individual_parameters_from_fit(fit, parameters, nonfixed, as_eta)
  
  # Check that parameters are correctly transformed with negative coefficients
  expect_equal(result$CL, 10 * exp(-0.1))
  expect_equal(result$V, 20 * exp(-0.2))
})

test_that("get_individual_parameters_from_fit handles zero coefficients", {
  # Create test data with zero coefficients
  fit <- list(coef = c(0, 0))
  parameters <- list(CL = 10, V = 20)
  nonfixed <- c("CL", "V")
  as_eta <- NULL
  
  result <- get_individual_parameters_from_fit(fit, parameters, nonfixed, as_eta)
  
  # Check that parameters remain unchanged with zero coefficients
  expect_equal(result$CL, 10)
  expect_equal(result$V, 20)
})

test_that("get_individual_parameters_from_fit handles mixed parameter types", {
  # Create test data with mixed parameter types (some as_eta, some transformed)
  fit <- list(coef = c(0.1, 0.2, 0.3))
  parameters <- list(CL = 10, V = 20, KA = 30)
  nonfixed <- c("CL", "V", "KA")
  as_eta <- c("V")  # Only V is direct eta
  
  result <- get_individual_parameters_from_fit(fit, parameters, nonfixed, as_eta)
  
  # Check that parameters are correctly handled
  expect_equal(result$CL, 10 * exp(0.1))
  expect_equal(result$V, 0.2)
  expect_equal(result$KA, 30 * exp(0.3))
}) 