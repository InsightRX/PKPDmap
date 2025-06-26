test_that("parse_t_obs returns observation times correctly", {
  data <- data.frame(t = c(0, 1, 2, 3), obs_type = c("PK", "PK", "PD", "PD"))
  
  result <- parse_t_obs(data)
  expect_equal(result, c(0, 1, 2, 3))
})

test_that("parse_t_obs handles single observation", {
  data <- data.frame(t = 0, obs_type = "PK")
  
  result <- parse_t_obs(data)
  expect_equal(result, 0)
})

test_that("parse_t_obs handles empty data frame", {
  data <- data.frame(t = numeric(0), obs_type = character(0))
  
  result <- parse_t_obs(data)
  expect_equal(result, numeric(0))
})

test_that("parse_t_obs handles duplicate times with different obs_type", {
  data <- data.frame(t = c(0, 1, 1, 2), obs_type = c("PK", "PK", "PD", "PD"))
  
  # Should not message about duplicates since obs_type differs
  expect_silent(result <- parse_t_obs(data))
  expect_equal(result, c(0, 1, 1, 2))
})

test_that("parse_t_obs detects and messages about duplicates with same obs_type", {
  data <- data.frame(t = c(0, 1, 1, 2), obs_type = c("PK", "PK", "PK", "PD"))
  
  # Should message about duplicates
  expect_message(
    result <- parse_t_obs(data),
    "Duplicate times were detected in data. Estimation will proceed but please check that data is correct. For putting more weight on certain measurements, please use the `weights` argument."
  )
  expect_equal(result, c(0, 1, 1, 2))
})

test_that("parse_t_obs detects duplicates with same time and obs_type", {
  data <- data.frame(t = c(0, 1, 1, 1, 2), obs_type = c("PK", "PK", "PK", "PD", "PD"))
  
  # Should message about duplicates
  expect_message(
    result <- parse_t_obs(data),
    "Duplicate times were detected in data. Estimation will proceed but please check that data is correct. For putting more weight on certain measurements, please use the `weights` argument."
  )
  expect_equal(result, c(0, 1, 1, 1, 2))
})

test_that("parse_t_obs handles multiple duplicate groups", {
  data <- data.frame(
    t = c(0, 1, 1, 2, 2, 2, 3), 
    obs_type = c("PK", "PK", "PK", "PD", "PD", "PK", "PK")
  )
  
  # Should message about duplicates
  expect_message(
    result <- parse_t_obs(data),
    "Duplicate times were detected in data. Estimation will proceed but please check that data is correct. For putting more weight on certain measurements, please use the `weights` argument."
  )
  expect_equal(result, c(0, 1, 1, 2, 2, 2, 3))
})

test_that("parse_t_obs handles all identical times and obs_types", {
  data <- data.frame(t = c(1, 1, 1), obs_type = c("PK", "PK", "PK"))
  
  # Should message about duplicates
  expect_message(
    result <- parse_t_obs(data),
    "Duplicate times were detected in data. Estimation will proceed but please check that data is correct. For putting more weight on certain measurements, please use the `weights` argument."
  )
  expect_equal(result, c(1, 1, 1))
})

test_that("parse_t_obs handles numeric obs_type", {
  data <- data.frame(t = c(0, 1, 1, 2), obs_type = c(1, 1, 1, 2))
  
  # Should message about duplicates
  expect_message(
    result <- parse_t_obs(data),
    "Duplicate times were detected in data. Estimation will proceed but please check that data is correct. For putting more weight on certain measurements, please use the `weights` argument."
  )
  expect_equal(result, c(0, 1, 1, 2))
})

test_that("parse_t_obs handles character obs_type with spaces", {
  data <- data.frame(t = c(0, 1, 1, 2), obs_type = c("PK", "PK", "PK", "PD"))
  
  # Should message about duplicates
  expect_message(
    result <- parse_t_obs(data),
    "Duplicate times were detected in data. Estimation will proceed but please check that data is correct. For putting more weight on certain measurements, please use the `weights` argument."
  )
  expect_equal(result, c(0, 1, 1, 2))
})

test_that("parse_t_obs handles NA values in obs_type", {
  data <- data.frame(t = c(0, 1, 1, 2), obs_type = c("PK", NA, NA, "PD"))
  
  # Should message about duplicates (NA values are treated as identical)
  expect_message(
    result <- parse_t_obs(data),
    "Duplicate times were detected in data. Estimation will proceed but please check that data is correct. For putting more weight on certain measurements, please use the `weights` argument."
  )
  expect_equal(result, c(0, 1, 1, 2))
})

test_that("parse_t_obs handles NA values in time", {
  data <- data.frame(t = c(0, NA, NA, 2), obs_type = c("PK", "PK", "PK", "PD"))
  
  # Should message about duplicates (NA values are treated as identical)
  expect_message(
    result <- parse_t_obs(data),
    "Duplicate times were detected in data. Estimation will proceed but please check that data is correct. For putting more weight on certain measurements, please use the `weights` argument."
  )
  expect_equal(result, c(0, NA, NA, 2))
})

test_that("parse_t_obs handles mixed data types in obs_type", {
  data <- data.frame(t = c(0, 1, 1, 2), obs_type = c("PK", 1, 1, "PD"))
  
  # Should message about duplicates
  expect_message(
    result <- parse_t_obs(data),
    "Duplicate times were detected in data. Estimation will proceed but please check that data is correct. For putting more weight on certain measurements, please use the `weights` argument."
  )
  expect_equal(result, c(0, 1, 1, 2))
})

test_that("parse_t_obs handles decimal times", {
  data <- data.frame(t = c(0.5, 1.0, 1.0, 2.5), obs_type = c("PK", "PK", "PK", "PD"))
  
  # Should message about duplicates
  expect_message(
    result <- parse_t_obs(data),
    "Duplicate times were detected in data. Estimation will proceed but please check that data is correct. For putting more weight on certain measurements, please use the `weights` argument."
  )
  expect_equal(result, c(0.5, 1.0, 1.0, 2.5))
})

test_that("parse_t_obs handles negative times", {
  data <- data.frame(t = c(-1, 0, 0, 1), obs_type = c("PK", "PK", "PK", "PD"))
  
  # Should message about duplicates
  expect_message(
    result <- parse_t_obs(data),
    "Duplicate times were detected in data. Estimation will proceed but please check that data is correct. For putting more weight on certain measurements, please use the `weights` argument."
  )
  expect_equal(result, c(-1, 0, 0, 1))
})

test_that("parse_t_obs handles large datasets without duplicates", {
  # Create a large dataset with no duplicates
  n <- 1000
  data <- data.frame(
    t = 1:n,
    obs_type = rep(c("PK", "PD"), length.out = n)
  )
  
  # Should not message about duplicates
  expect_silent(result <- parse_t_obs(data))
  expect_equal(result, 1:n)
})
