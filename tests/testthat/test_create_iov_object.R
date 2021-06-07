
test_that("when no IOV present should throw error", {
 expect_error(
   create_iov_object(
     cv = list(CL = 0.201),
     omega = c(0.177241, 0, 0.992016, 0, 0, 0.725904, 0, 0, 0, 1),
     bins = c(0L, 24L, 48L, 72L, 9999L),
     parameters = structure(
       list(CL = 36.6, V = 496, Q = 31.7, Vp = 1270, KA = 1.7, TLAG = 0.41), 
       units = ""
     ),
     fixed = c("Q", "TLAG"),
     ruv = list(prop = 0.23, add = 0.87)
   )
 )
})

test_that("non-ordered fixed parameters are sorted properly", {
  res1 <- create_iov_object(
    cv = list(CL = 0.201),
    omega = c(0.177241, 0, 0.992016, 0, 0, 0.725904, 0, 0, 0, 1),
    bins = c(0L, 24L, 48L, 72L, 9999L),
    parameters = structure(
      list(
        CL = 36.6, V = 496, Q = 31.7, Vp = 1270, KA = 1.7, TLAG = 0.41, 
        kappa_CL_1 = 0, kappa_CL_2 = 0, kappa_CL_3 = 0, kappa_CL_4 = 0
      ), 
      units = ""
    ),
    fixed = c("Q", "TLAG"),
    ruv = list(prop = 0.23, add = 0.87)
  )
  expect_equal(
    names(res1$parameters), 
    c(
      "CL", "V", "Vp", "KA", 
      "kappa_CL_1", "kappa_CL_2", "kappa_CL_3", "kappa_CL_4", 
      "Q", "TLAG"
    )
  )
  expect_equal(length(res1$omega), 36)
})

test_that("With IOV but no fixed par, assume all in correct order", {
  res2 <- create_iov_object(
    cv = list(CL = 0.201),
    omega = c(0.177241, 0, 0.992016, 0, 0, 0.725904, 0, 0, 0, 1),
    bins = c(0L, 24L, 48L, 72L, 9999L),
    parameters = structure(
      list(
        CL = 36.6, V = 496, Q = 31.7, Vp = 1270, KA = 1.7, TLAG = 0.41, 
        kappa_CL_1 = 0, kappa_CL_2 = 0, kappa_CL_3 = 0, kappa_CL_4 = 0
      ), 
      units = ""
    ),
    ruv = list(prop = 0.23, add = 0.87)
  )
  expect_equal(
    names(res2$parameters), 
    c(
      "CL", "V", "Q", "Vp", 
      "kappa_CL_1", "kappa_CL_2", "kappa_CL_3", "kappa_CL_4", 
      "KA", "TLAG"
    )
  )
  expect_equal(length(res2$omega), 36)
})


test_that("Extraneous fixed parameters are ignored", {
  res1 <- create_iov_object(
    cv = list(CL = 0.201),
    omega = c(0.177241, 0, 0.992016),
    bins = c(0L, 24L, 48L, 72L, 9999L),
    parameters = structure(
      list(CL = 36.6, V = 496, TLAG = 0.41, kappa_CL_1 = 0, kappa_CL_2 = 0, kappa_CL_3 = 0, kappa_CL_4 = 0), 
      units = ""
    ),
    fixed = c("TLAG", "ForeignParameter"),
    ruv = list(prop = 0.23, add = 0.87)
  )
  expect_equal(
    names(res1$parameters), 
    c(
      "CL", "V",
      "kappa_CL_1", "kappa_CL_2", "kappa_CL_3", "kappa_CL_4", 
      "TLAG"
    )
  )
})

