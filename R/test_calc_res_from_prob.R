test_that("Calculation of residuals works when only bloq, both single and multiple", {
  p <- seq(0.3)
  expect_equal(
    round(calc_res_from_prob(p), 3),
    0.045
  )
  p <- seq(0, 1, .1)
  expect_equal(
    round(calc_res_from_prob(p), 3),
    c(4.292, 2.146, 1.794, 1.552, 1.354, 1.177, 1.011, 0.845, 0.668, 
      0.459, 0.045)
  )
})

test_that("Calculation of residuals doesn't error on NA or NULL", {
  expect_equal(
    calc_res_from_prob(NA),
    NA_real_
  )
  expect_equal(
    calc_res_from_prob(NaN),
    NA_real_
  )
  expect_equal(
    length(calc_res_from_prob(NULL)),
    0
  )
})

