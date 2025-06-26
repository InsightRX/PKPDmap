test_that("parse_fixed returns correct fixed parameters when all exist", {
  parameters <- list(CL = 1, V = 10, KA = 2)
  fixed <- c("CL", "V")
  
  result <- parse_fixed(fixed, parameters)
  expect_equal(result, c("CL", "V"))
})

test_that("parse_fixed returns correct fixed parameters when some exist", {
  parameters <- list(CL = 1, V = 10, KA = 2)
  fixed <- c("CL", "V", "NONEXISTENT")
  
  # Should warn about missing parameters
  expect_warning(
    result <- parse_fixed(fixed, parameters),
    "Warning: not all fixed parameters were found in parameter set!"
  )
  
  # Should return only existing parameters
  expect_equal(result, c("CL", "V"))
})

test_that("parse_fixed returns NULL when no fixed parameters exist", {
  parameters <- list(CL = 1, V = 10, KA = 2)
  fixed <- c("NONEXISTENT1", "NONEXISTENT2")
  
  # Should warn about missing parameters
  expect_warning(
    result <- parse_fixed(fixed, parameters),
    "Warning: not all fixed parameters were found in parameter set!"
  )
  
  # Should return NULL when no parameters match
  expect_null(result)
})

test_that("parse_fixed returns NULL when fixed is empty", {
  parameters <- list(CL = 1, V = 10, KA = 2)
  fixed <- character(0)
  
  result <- parse_fixed(fixed, parameters)
  expect_null(result)
})

test_that("parse_fixed handles case sensitivity correctly", {
  parameters <- list(CL = 1, V = 10, KA = 2)
  fixed <- c("cl", "V")  # Note: "cl" is lowercase
  
  # Should warn about missing parameters (case sensitive)
  expect_warning(
    result <- parse_fixed(fixed, parameters),
    "Warning: not all fixed parameters were found in parameter set!"
  )
  
  # Should return only exact matches
  expect_equal(result, "V")
})

test_that("parse_fixed handles single parameter correctly", {
  parameters <- list(CL = 1, V = 10, KA = 2)
  fixed <- "CL"
  
  result <- parse_fixed(fixed, parameters)
  expect_equal(result, "CL")
})

test_that("parse_fixed handles single non-existent parameter", {
  parameters <- list(CL = 1, V = 10, KA = 2)
  fixed <- "NONEXISTENT"
  
  # Should warn about missing parameter
  expect_warning(
    result <- parse_fixed(fixed, parameters),
    "Warning: not all fixed parameters were found in parameter set!"
  )
  
  # Should return NULL
  expect_null(result)
})

test_that("parse_fixed handles parameters with special characters", {
  parameters <- list("CL_1" = 1, "V_central" = 10, "KA_abs" = 2)
  fixed <- c("CL_1", "V_central")
  
  result <- parse_fixed(fixed, parameters)
  expect_equal(result, c("CL_1", "V_central"))
})

test_that("parse_fixed handles mixed existing and non-existing parameters", {
  parameters <- list(CL = 1, V = 10, KA = 2)
  fixed <- c("CL", "NONEXISTENT", "V")
  
  # Should warn about missing parameters
  expect_warning(
    result <- parse_fixed(fixed, parameters),
    "Warning: not all fixed parameters were found in parameter set!"
  )
  
  # Should return existing parameters in order they appear in fixed
  expect_equal(result, c("CL", "V"))
})

test_that("parse_fixed handles duplicate parameters in fixed", {
  parameters <- list(CL = 1, V = 10, KA = 2)
  fixed <- c("CL", "V", "CL")  # CL appears twice
  
  result <- parse_fixed(fixed, parameters)
  expect_equal(result, c("CL", "V"))
})

test_that("parse_fixed handles empty parameters list", {
  parameters <- list()
  fixed <- c("CL", "V")
  
  # Should warn about missing parameters
  expect_warning(
    result <- parse_fixed(fixed, parameters),
    "Warning: not all fixed parameters were found in parameter set!"
  )
  
  # Should return NULL
  expect_null(result)
})

test_that("parse_fixed handles NULL fixed parameter", {
  parameters <- list(CL = 1, V = 10, KA = 2)
  
  # NULL should be treated as empty character vector
  expect_null(parse_fixed(NULL, parameters))
})

test_that("parse_fixed handles numeric parameters", {
  parameters <- list(CL = 1, V = 10, KA = 2)
  fixed <- c("CL", "V")
  
  result <- parse_fixed(fixed, parameters)
  expect_equal(result, c("CL", "V"))
})

test_that("parse_fixed handles character parameters", {
  parameters <- list(CL = "value1", V = "value2", KA = "value3")
  fixed <- c("CL", "V")
  
  result <- parse_fixed(fixed, parameters)
  expect_equal(result, c("CL", "V"))
}) 

