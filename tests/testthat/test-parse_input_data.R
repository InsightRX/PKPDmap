test_that("parse_input_data handles regular data frame input", {
  # Create test data
  data <- data.frame(
    T = c(2, 1, 3),
    OBS_TYPE = c(2, 1, 1),
    Y = c(10, 20, 30)
  )
  
  result <- parse_input_data(data)
  
  # Check column names are lowercase
  expect_equal(names(result), c("t", "obs_type", "y", "obs_type"))
  
  # Check sorting
  expect_equal(result$t, c(1, 2, 3))
  expect_equal(result$obs_type, c(1, 2, 1))
})

test_that("parse_input_data handles PKPDsim object", {
  # Create mock PKPDsim object
  data <- structure(
    data.frame(
      t = c(2, 1, 3),
      comp = c("obs", "dose", "obs"),
      dv = c(10, 20, 30)
    ),
    class = c("data.frame", "PKPDsim")
  )
  
  result <- parse_input_data(data, cols=list(x = "t", y = "dv"))
  
  # Check that only obs rows are kept
  expect_equal(nrow(result), 2)
  expect_true(all(result$comp == "obs"))
  
  # Check that evid column is added and set to 0
  expect_true("evid" %in% names(result))
  expect_true(all(result$evid == 0))
})

test_that("parse_input_data handles EVID filtering", {
  # Create test data with EVID column
  data <- data.frame(
    t = c(2, 1, 3),
    obs_type = c(2, 1, 1),
    dv = c(10, 20, 30),
    EVID = c(0, 1, 0)
  )
  
  result <- parse_input_data(data, cols=list(x = "t", y = "dv"))
  
  # Check that only EVID=0 rows are kept
  expect_equal(nrow(result), 2)
  expect_true(all(result$evid == 0))
})

test_that("parse_input_data handles obs_type_label argument", {
  # Test with NULL obs_type_label (default behavior)
  data <- data.frame(
    t = c(2, 1, 3),
    dv = c(10, 20, 30)
  )
  result <- parse_input_data(data, cols=list(x = "t", y = "dv"))
  expect_equal(result$obs_type, c(1, 1, 1))
  
  # Test with custom obs_type_label
  data <- data.frame(
    t = c(2, 1, 3),
    dv = c(10, 20, 30),
    measurement_type = c(2, 1, 3)
  )
  result <- parse_input_data(
    data, 
    cols=list(x = "t", y = "dv"),
    obs_type_label = "measurement_type"
  )
  expect_equal(result$obs_type, c(1, 2, 3))
})

test_that("parse_input_data handles PKPDsim object with obs_type_label", {
  # Create mock PKPDsim object with obs_type_label
  data <- structure(
    data.frame(
      t = c(2, 1, 3),
      comp = c("obs", 1, "obs"),
      dv = c(10, 20, 30),
      measurement_type = c(2, 1, 3)
    ),
    class = c("data.frame", "PKPDsim")
  )
  
  result <- parse_input_data(
    data, 
    cols=list(x = "t", y = "dv"), 
    obs_type_label = "measurement_type"
  )
  
  # Check that only obs rows are kept and obs_type is set correctly
  expect_equal(nrow(result), 2)
  expect_true(all(result$comp == "obs"))
  expect_equal(result$obs_type, c(2, 3))
})

test_that("parse_input_data handles multiple EVID values", {
  # Create test data with multiple EVID values
  data <- data.frame(
    t = c(2, 1, 3, 4),
    obs_type = c(2, 1, 1, 1),
    dv = c(10, 20, 30, 40),
    EVID = c(0, 1, 2, 0)
  )
  
  result <- parse_input_data(data, cols=list(x = "t", y = "dv"))
  
  # Check that only EVID=0 rows are kept
  expect_equal(nrow(result), 2)
  expect_true(all(result$evid == 0))
  expect_equal(result$t, c(2, 4))
})

test_that("parse_input_data preserves data values after transformations", {
  # Create test data
  data <- data.frame(
    t = c(2, 1, 3),
    obs_type = c(2, 1, 1),
    DV = c(10, 20, 30),
    EVID = c(0, 1, 0)
  )
  
  result <- parse_input_data(data, cols=list(x = "t", y = "dv"))
  
  # Check that data values are preserved after transformations
  expect_equal(result$dv, c(10, 30))
  expect_equal(result$t, c(2, 3))
})

