test_that("parse_weight_prior handles default values", {
  # Test default weight_prior
  result <- parse_weight_prior()
  expect_equal(result, 1)  # 1^2 = 1
  
  # Test default weight_prior with explicit type
  result <- parse_weight_prior(type = "map")
  expect_equal(result, 1)
})

test_that("parse_weight_prior handles NULL and NA inputs", {
  # Test NULL weight_prior
  result <- parse_weight_prior(NULL)
  expect_equal(result, 1)
  
  # Test NA weight_prior
  result <- parse_weight_prior(NA)
  expect_equal(result, 1)
  
  # Test NULL weight_prior with explicit type
  result <- parse_weight_prior(NULL, type = "map")
  expect_equal(result, 1)
})

test_that("parse_weight_prior handles different weight values", {
  # Test with weight = 2
  result <- parse_weight_prior(2)
  expect_equal(result, 4)  # 2^2 = 4
  
  # Test with weight = 0.5
  result <- parse_weight_prior(0.5)
  expect_equal(result, 0.25)  # 0.5^2 = 0.25
  
  # Test with weight = 0
  result <- parse_weight_prior(0)
  expect_equal(result, 0)  # 0^2 = 0
})

test_that("parse_weight_prior handles PLS type", {
  # Test with lowercase "pls"
  result <- parse_weight_prior(type = "pls")
  expect_equal(result, 0.001)
  
  # Test with uppercase "PLS"
  result <- parse_weight_prior(type = "PLS")
  expect_equal(result, 0.001)
  
  # Test with mixed case "Pls"
  result <- parse_weight_prior(type = "Pls")
  expect_equal(result, 0.001)
  
  # Test with weight_prior and PLS type
  result <- parse_weight_prior(2, type = "pls")
  expect_equal(result, 0.001)  # Should ignore weight_prior for PLS
})

test_that("parse_weight_prior handles MAP type with different cases", {
  # Test with lowercase "map"
  result <- parse_weight_prior(2, type = "map")
  expect_equal(result, 4)
  
  # Test with uppercase "MAP"
  result <- parse_weight_prior(2, type = "MAP")
  expect_equal(result, 4)
  
  # Test with mixed case "Map"
  result <- parse_weight_prior(2, type = "Map")
  expect_equal(result, 4)
}) 