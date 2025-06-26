test_that("parse_error handles NULL input", {
  expect_error(
    parse_error(NULL),
    "No residual error specified, cannot perform fit."
  )
})

test_that("parse_error handles missing error components", {
  # Test with only prop specified
  error <- list(prop = 0.1)
  result <- parse_error(error)
  expect_equal(result$prop, 0.1)
  expect_equal(result$add, 0)
  expect_equal(result$exp, 0)
  
  # Test with only add specified
  error <- list(add = 0.2)
  result <- parse_error(error)
  expect_equal(result$prop, 0)
  expect_equal(result$add, 0.2)
  expect_equal(result$exp, 0)
  
  # Test with only exp specified
  error <- list(exp = 0.3)
  result <- parse_error(error)
  expect_equal(result$prop, 0)
  expect_equal(result$add, 0)
  expect_equal(result$exp, 0.3)
})

test_that("parse_error handles all zeros", {
  error <- list(prop = 0, add = 0, exp = 0)
  expect_error(
    parse_error(error),
    "No residual error model specified, or residual error is 0."
  )
})

test_that("parse_error handles valid error specifications", {
  # Test with all components specified
  error <- list(prop = 0.1, add = 0.2, exp = 0.3)
  result <- parse_error(error)
  expect_equal(result$prop, 0.1)
  expect_equal(result$add, 0.2)
  expect_equal(result$exp, 0.3)
  
  # Test with mixed zero and non-zero components
  error <- list(prop = 0.1, add = 0, exp = 0.3)
  result <- parse_error(error)
  expect_equal(result$prop, 0.1)
  expect_equal(result$add, 0)
  expect_equal(result$exp, 0.3)
}) 