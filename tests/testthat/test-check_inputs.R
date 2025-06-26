test_that("check_inputs passes with valid MAP inputs", {
  # Create mock objects
  model <- function() {}
  data <- data.frame(time = 1:3, dv = c(1, 2, 3))
  parameters <- list(CL = 1, V = 10)
  omega <- matrix(c(0.1, 0, 0, 0.1), nrow = 2)
  regimen <- list(dose = 100, interval = 12)
  
  # Should not throw any errors
  expect_no_error(check_inputs(model, data, parameters, omega, regimen, NULL, "MAP"))
})

test_that("check_inputs passes with valid PLS inputs", {
  # Create mock objects
  model <- function() {}
  data <- data.frame(time = 1:3, dv = c(1, 2, 3))
  parameters <- list(CL = 1, V = 10)
  omega <- matrix(c(0.1, 0, 0, 0.1), nrow = 2)
  regimen <- list(dose = 100, interval = 12)
  
  # Should not throw any errors
  expect_no_error(check_inputs(model, data, parameters, omega, regimen, NULL, "PLS"))
})

test_that("check_inputs passes with valid pls inputs (lowercase)", {
  # Create mock objects
  model <- function() {}
  data <- data.frame(time = 1:3, dv = c(1, 2, 3))
  parameters <- list(CL = 1, V = 10)
  omega <- matrix(c(0.1, 0, 0, 0.1), nrow = 2)
  regimen <- list(dose = 100, interval = 12)
  
  # Should not throw any errors
  expect_no_error(check_inputs(model, data, parameters, omega, regimen, NULL, "pls"))
})

test_that("check_inputs passes with valid MAP inputs (lowercase)", {
  # Create mock objects
  model <- function() {}
  data <- data.frame(time = 1:3, dv = c(1, 2, 3))
  parameters <- list(CL = 1, V = 10)
  omega <- matrix(c(0.1, 0, 0, 0.1), nrow = 2)
  regimen <- list(dose = 100, interval = 12)
  
  # Should not throw any errors
  expect_no_error(check_inputs(model, data, parameters, omega, regimen, NULL, "map"))
})

test_that("check_inputs passes with valid censoring argument", {
  # Create mock objects
  model <- function() {}
  data <- data.frame(time = 1:3, dv = c(1, 2, 3))
  parameters <- list(CL = 1, V = 10)
  omega <- matrix(c(0.1, 0, 0, 0.1), nrow = 2)
  regimen <- list(dose = 100, interval = 12)
  
  # Should not throw any errors with character censoring
  expect_no_error(check_inputs(model, data, parameters, omega, regimen, "cens_col", "MAP"))
})

test_that("check_inputs passes with other type values", {
  # Create mock objects
  model <- function() {}
  data <- data.frame(time = 1:3, dv = c(1, 2, 3))
  parameters <- list(CL = 1, V = 10)
  omega <- matrix(c(0.1, 0, 0, 0.1), nrow = 2)
  regimen <- list(dose = 100, interval = 12)
  
  # Should not throw any errors for non-MAP/PLS types
  expect_no_error(check_inputs(model, data, parameters, omega, regimen, NULL, "other_type"))
})

test_that("check_inputs fails when model is NULL for MAP type", {
  data <- data.frame(time = 1:3, dv = c(1, 2, 3))
  parameters <- list(CL = 1, V = 10)
  omega <- matrix(c(0.1, 0, 0, 0.1), nrow = 2)
  regimen <- list(dose = 100, interval = 12)
  
  expect_error(
    check_inputs(NULL, data, parameters, omega, regimen, NULL, "MAP"),
    "The 'model', 'data', 'omega', 'regimen', and 'parameters' arguments are required."
  )
})

test_that("check_inputs fails when data is NULL for MAP type", {
  model <- function() {}
  parameters <- list(CL = 1, V = 10)
  omega <- matrix(c(0.1, 0, 0, 0.1), nrow = 2)
  regimen <- list(dose = 100, interval = 12)
  
  expect_error(
    check_inputs(model, NULL, parameters, omega, regimen, NULL, "MAP"),
    "The 'model', 'data', 'omega', 'regimen', and 'parameters' arguments are required."
  )
})

test_that("check_inputs fails when parameters is NULL for MAP type", {
  model <- function() {}
  data <- data.frame(time = 1:3, dv = c(1, 2, 3))
  omega <- matrix(c(0.1, 0, 0, 0.1), nrow = 2)
  regimen <- list(dose = 100, interval = 12)
  
  expect_error(
    check_inputs(model, data, NULL, omega, regimen, NULL, "MAP"),
    "The 'model', 'data', 'omega', 'regimen', and 'parameters' arguments are required."
  )
})

test_that("check_inputs fails when omega is NULL for MAP type", {
  model <- function() {}
  data <- data.frame(time = 1:3, dv = c(1, 2, 3))
  parameters <- list(CL = 1, V = 10)
  regimen <- list(dose = 100, interval = 12)
  
  expect_error(
    check_inputs(model, data, parameters, NULL, regimen, NULL, "MAP"),
    "The 'model', 'data', 'omega', 'regimen', and 'parameters' arguments are required."
  )
})

test_that("check_inputs fails when regimen is NULL for MAP type", {
  model <- function() {}
  data <- data.frame(time = 1:3, dv = c(1, 2, 3))
  parameters <- list(CL = 1, V = 10)
  omega <- matrix(c(0.1, 0, 0, 0.1), nrow = 2)
  
  expect_error(
    check_inputs(model, data, parameters, omega, NULL, NULL, "MAP"),
    "The 'model', 'data', 'omega', 'regimen', and 'parameters' arguments are required."
  )
})

test_that("check_inputs fails when model is not a function", {
  model <- "not_a_function"
  data <- data.frame(time = 1:3, dv = c(1, 2, 3))
  parameters <- list(CL = 1, V = 10)
  omega <- matrix(c(0.1, 0, 0, 0.1), nrow = 2)
  regimen <- list(dose = 100, interval = 12)
  
  expect_error(
    check_inputs(model, data, parameters, omega, regimen, NULL, "MAP"),
    "The 'model' argument requires a function, e.g. a model defined using the new_ode_model\\(\\) function from the PKPDsim package."
  )
})

test_that("check_inputs fails when model is a list (not a function)", {
  model <- list(not_a_function = TRUE)
  data <- data.frame(time = 1:3, dv = c(1, 2, 3))
  parameters <- list(CL = 1, V = 10)
  omega <- matrix(c(0.1, 0, 0, 0.1), nrow = 2)
  regimen <- list(dose = 100, interval = 12)
  
  expect_error(
    check_inputs(model, data, parameters, omega, regimen, NULL, "MAP"),
    "The 'model' argument requires a function, e.g. a model defined using the new_ode_model\\(\\) function from the PKPDsim package."
  )
})

test_that("check_inputs fails when censoring is not NULL and not character", {
  model <- function() {}
  data <- data.frame(time = 1:3, dv = c(1, 2, 3))
  parameters <- list(CL = 1, V = 10)
  omega <- matrix(c(0.1, 0, 0, 0.1), nrow = 2)
  regimen <- list(dose = 100, interval = 12)
  
  # Test with numeric censoring
  expect_error(
    check_inputs(model, data, parameters, omega, regimen, 123, "MAP"),
    "Censoring argument requires label specifying column in dataset with censoring info."
  )
  
  # Test with logical censoring
  expect_error(
    check_inputs(model, data, parameters, omega, regimen, TRUE, "MAP"),
    "Censoring argument requires label specifying column in dataset with censoring info."
  )
  
  # Test with list censoring
  expect_error(
    check_inputs(model, data, parameters, omega, regimen, list(cens = "col"), "MAP"),
    "Censoring argument requires label specifying column in dataset with censoring info."
  )
})

test_that("check_inputs allows NULL censoring", {
  model <- function() {}
  data <- data.frame(time = 1:3, dv = c(1, 2, 3))
  parameters <- list(CL = 1, V = 10)
  omega <- matrix(c(0.1, 0, 0, 0.1), nrow = 2)
  regimen <- list(dose = 100, interval = 12)
  
  # Should not throw any errors with NULL censoring
  expect_no_error(check_inputs(model, data, parameters, omega, regimen, NULL, "MAP"))
})

test_that("check_inputs allows character censoring", {
  model <- function() {}
  data <- data.frame(time = 1:3, dv = c(1, 2, 3))
  parameters <- list(CL = 1, V = 10)
  omega <- matrix(c(0.1, 0, 0, 0.1), nrow = 2)
  regimen <- list(dose = 100, interval = 12)
  
  # Should not throw any errors with character censoring
  expect_no_error(check_inputs(model, data, parameters, omega, regimen, "censoring_column", "MAP"))
})

test_that("check_inputs works with different function types", {
  # Test with anonymous function
  model <- function(x) x + 1
  data <- data.frame(time = 1:3, dv = c(1, 2, 3))
  parameters <- list(CL = 1, V = 10)
  omega <- matrix(c(0.1, 0, 0, 0.1), nrow = 2)
  regimen <- list(dose = 100, interval = 12)
  
  expect_no_error(check_inputs(model, data, parameters, omega, regimen, NULL, "MAP"))
  
  # Test with named function
  test_model <- function() {}
  expect_no_error(check_inputs(test_model, data, parameters, omega, regimen, NULL, "MAP"))
}) 