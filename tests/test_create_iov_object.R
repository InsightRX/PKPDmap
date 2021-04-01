library(PKPDmap)
library(testit)

## Test output from create_iov_object, when no IOV present (should throw error)
assert("Throw error when no IOV present", has_error(create_iov_object(
  cv = list(CL = 0.201),
  omega = c(0.177241, 0, 0.992016, 0, 0, 0.725904, 0, 0, 0, 1),
  bins = c(0L, 24L, 48L, 72L, 9999L),
  parameters = structure(list(CL = 36.6, V = 496, Q = 31.7, Vp = 1270, KA = 1.7, TLAG = 0.41), units = ""),
  fixed = c("Q", "TLAG"),
  ruv = list(prop = 0.23, add = 0.87)
), silent = TRUE))

## Test output from create_iov_object, with non-ordered fixed parameters
res1 <- create_iov_object(
  cv = list(CL = 0.201),
  omega = c(0.177241, 0, 0.992016, 0, 0, 0.725904, 0, 0, 0, 1),
  bins = c(0L, 24L, 48L, 72L, 9999L),
  parameters = structure(list(CL = 36.6, V = 496, Q = 31.7, Vp = 1270, KA = 1.7, TLAG = 0.41, kappa_CL_1 = 0, kappa_CL_2 = 0, kappa_CL_3 = 0, kappa_CL_4 = 0), units = ""),
  fixed = c("Q", "TLAG"),
  ruv = list(prop = 0.23, add = 0.87)
)

assert("Proper order of parameters to match omega", all(names(res1$parameters) == c("CL", "V", "Vp", "KA", "kappa_CL_1", "kappa_CL_2", "kappa_CL_3", "kappa_CL_4", "Q", "TLAG")))
assert("Correct length of omega", length(res1$omega) == 36)

## Test output from create_iov_object, with IOV but no specified fixed parameters
## in this case it will just assume the user specified parameters in correct order.
res2 <- create_iov_object(
  cv = list(CL = 0.201),
  omega = c(0.177241, 0, 0.992016, 0, 0, 0.725904, 0, 0, 0, 1),
  bins = c(0L, 24L, 48L, 72L, 9999L),
  parameters = structure(list(CL = 36.6, V = 496, Q = 31.7, Vp = 1270, KA = 1.7, TLAG = 0.41, kappa_CL_1 = 0, kappa_CL_2 = 0, kappa_CL_3 = 0, kappa_CL_4 = 0), units = ""),
  ruv = list(prop = 0.23, add = 0.87)
)

assert("Proper order of parameters to match omega", all(names(res2$parameters) == c("CL", "V", "Q", "Vp", "kappa_CL_1", "kappa_CL_2", "kappa_CL_3", "kappa_CL_4", "KA", "TLAG")))
assert("Correct length of omega", length(res1$omega) == 36)
