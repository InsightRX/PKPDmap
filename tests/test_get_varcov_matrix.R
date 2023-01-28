fallback <- matrix(c(
  0.2, 0.02,
  0.02, 0.2
), ncol = 2)

test_that("get_varcov_matrix returns matrix if present and OK", {
  expect_equal(
    get_varcov_matrix(list(
      vcov_full = matrix(c(
        0.1, 0.01,
        0.01, 0.1
      ), ncol = 2)
    ), fallback = fallback),
    matrix(
      c(0.1, 0.01, 0.01, 0.1),
      ncol = 2
    )
  )
})

test_that("get_varcov_matrix returns fallback if not present", {
  expect_equal(
    get_varcov_matrix(NULL, fallback = fallback),
    fallback
  )
})

test_that("get_varcov_matrix returns fallback if not pos-def", {
  expect_warning(
    expect_equal(
      get_varcov_matrix(list(
        vcov_full = matrix(c(
          0.1, 0.2,
          0.2, 0.1
        ), ncol = 2)
      ), fallback = fallback),
      fallback
    )
  )
})

test_that("get_varcov_matrix returns fallback if has NAs", {
  expect_warning(
    expect_equal(
      get_varcov_matrix(list(
        vcov_full = matrix(c(
          0.1, NA,
          NA, 0.1
        ), ncol = 2)
      ), fallback = fallback),
      fallback
    )
  )
})

