test_that("parse_censoring handles NULL censoring", {
  data <- data.frame(dv = c(1, 2, 3))
  result <- parse_censoring(NULL, data)
  expect_null(result)
})

test_that("parse_censoring handles no censored values", {
  data <- data.frame(dv = c(1, 2, 3), cens = c(0, 0, 0))
  result <- parse_censoring("cens", data)
  expect_null(result)
  
  # Test with verbose messaging
  expect_message(
    parse_censoring("cens", data, verbose = TRUE),
    "Warning: censoring specified, but no censored values in data."
  )
})

test_that("parse_censoring handles censored values", {
  data <- data.frame(dv = c(1, 2, 3), cens = c(0, 1, 0))
  result <- parse_censoring("cens", data)
  expect_equal(result, c(FALSE, TRUE, FALSE))
  
  # Test with verbose messaging
  expect_message(
    parse_censoring("cens", data, verbose = TRUE),
    "One or more values in data are censored, including censoring in likelihood."
  )
})

test_that("parse_censoring is case insensitive", {  
  # Test with lowercase column name
  data <- data.frame(dv = c(1, 2, 3), cens = c(0, 1, 0))
  result <- parse_censoring("CENS", data)
  expect_equal(result, c(FALSE, TRUE, FALSE))
}) 
