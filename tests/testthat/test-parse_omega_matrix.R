test_that("parse_omega_matrix handles matrix input", {
  # Create test matrix
  omega <- matrix(c(0.1, 0.2, 0.2, 0.3), nrow = 2, ncol = 2)
  parameters <- c("CL", "V")
  fixed <- character(0)
  
  result <- parse_omega_matrix(omega, parameters, fixed)
  
  # Check that matrix is returned unchanged
  expect_equal(result, omega)
  expect_equal(dim(result), c(2, 2))
})

test_that("parse_omega_matrix handles triangle matrix input", {
  # Create test triangle matrix (lower triangle)
  omega <- c(0.1, 0.2, 0.3)  # represents matrix c(0.1, 0.2, 0.2, 0.3)
  parameters <- c("CL", "V")
  fixed <- character(0)
  
  result <- parse_omega_matrix(omega, parameters, fixed)
  
  # Check that triangle matrix is converted to full matrix
  expect_equal(dim(result), c(2, 2))
  expect_equal(result[1,1], 0.1)
  expect_equal(result[2,1], 0.2)
  expect_equal(result[1,2], 0.2)
  expect_equal(result[2,2], 0.3)
})

test_that("parse_omega_matrix validates matrix size correctly", {
  # Test with matching sizes
  omega <- c(0.1, 0.2, 0.3)  # 2x2 matrix
  parameters <- c("CL", "V")
  fixed <- character(0)
  
  result <- parse_omega_matrix(omega, parameters, fixed)
  expect_equal(dim(result), c(2, 2))
  
  # Test with fixed parameters
  omega <- c(0.1, 0.2, 0.3)  # 2x2 matrix
  parameters <- c("CL", "V", "KA")
  fixed <- "KA"
  
  result <- parse_omega_matrix(omega, parameters, fixed)
  expect_equal(dim(result), c(2, 2))
})

test_that("parse_omega_matrix handles error cases", {
  # Test with too small omega matrix
  omega <- c(0.1)  # 1x1 matrix
  parameters <- c("CL", "V")
  fixed <- character(0)
  
  expect_error(
    parse_omega_matrix(omega, parameters, fixed),
    "Provided omega matrix is smaller than expected based on the number of model parameters"
  )
  
  # Test with too small omega matrix and fixed parameters
  omega <- c(0.1)  # 1x1 matrix
  parameters <- c("CL", "V", "KA")
  fixed <- "KA"
  
  expect_error(
    parse_omega_matrix(omega, parameters, fixed),
    "Provided omega matrix is smaller than expected based on the number of model parameters"
  )
})

test_that("parse_omega_matrix handles edge cases", {
  # Test with single parameter
  omega <- c(0.1)  # 1x1 matrix
  parameters <- c("CL")
  fixed <- character(0)
  
  result <- parse_omega_matrix(omega, parameters, fixed)
  expect_equal(dim(result), c(1, 1))
  expect_equal(result[1,1], 0.1)
  
  # Test with all parameters fixed
  omega <- c(0.1)  # 1x1 matrix
  parameters <- c("CL", "V")
  fixed <- c("CL", "V")
  
  result <- parse_omega_matrix(omega, parameters, fixed)
  expect_equal(dim(result), c(1, 1))
}) 