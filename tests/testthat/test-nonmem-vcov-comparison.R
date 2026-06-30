# Cross-check the FOCE-style variance-covariance matrix produced by
# get_map_estimates() against NONMEM 7.5.1.
#
# The NONMEM reference values below are HARD-CODED: NONMEM is not run during
# CI (it is licensed and not installed on the CI image). The numbers were
# produced once with NONMEM 7.5.1 using the control stream documented below,
# and this test reruns only the PKPDmap side and compares against them.
#
# Both sides estimate per-subject etas with the population parameters fixed
# (NONMEM: $EST METHOD=COND INTERACTION MAXEVAL=0 POSTHOC), which is exactly
# what get_map_estimates() does (MAP / empirical Bayes with fixed thetas and
# omega). The FOCE conditional eta variance-covariance matrix NONMEM writes to
# its .phi file (the ETC(i,j) columns) is the same quantity PKPDmap returns in
# obj$vcov_full.
#
# INTERACTION is required on the NONMEM side: PKPDmap evaluates the residual
# error variance at the *individual* prediction (IPRED). Plain FOCE (no
# INTERACTION) evaluates it at the population prediction (eta = 0), which
# changes the eta mode for proportional-error models and would not match.
#
# Per-observation weights are checked with a non-uniform weight vector
# c(1, 3): the t = 8 observation is replicated three times in the NONMEM data
# set (and t = 2 left as a single record), an exact integer weighting of the
# likelihood contributions. Non-uniform weights exercise the per-observation
# diag(weights / sigma^2) term in the FOCE Hessian, not just a global scaling.
#
# --- NONMEM control stream (unit weights) -----------------------------------
# $PROBLEM 1cmt IV infusion, single subject, MAP (posthoc) vcov check
# $INPUT ID TIME AMT RATE DV MDV EVID CMT
# $DATA data.csv IGNORE=@
# $SUBROUTINES ADVAN1 TRANS2
# $PK
#   CL = THETA(1)*EXP(ETA(1))
#   V  = THETA(2)*EXP(ETA(2))
#   S1 = V
# $ERROR
#   IPRED = F
#   Y = IPRED*(1 + EPS(1)) + EPS(2)
# $THETA  5 FIX  50 FIX
# $OMEGA  0.1  0.1
# $SIGMA  0.01  0.01
# $ESTIMATION METHOD=COND INTERACTION MAXEVAL=0 POSTHOC NSIG=3 SIGL=9
# $COVARIANCE
#
# data.csv (unit weights):       weighted run, weights c(1, 3) (t=8 x3):
#   ID,TIME,AMT,RATE,DV,MDV,EVID,CMT      1,0,1000,1000,0,1,1,1
#   1,0,1000,1000,0,1,1,1                 1,2,0,0,18,0,0,1
#   1,2,0,0,18,0,0,1                      1,8,0,0,6,0,0,1
#   1,8,0,0,6,0,0,1                       1,8,0,0,6,0,0,1
#                                         1,8,0,0,6,0,0,1
# ----------------------------------------------------------------------------

test_that("FOCE vcov matches NONMEM 7.5.1 (unit and per-observation weights)", {
  # PKPDmap and NONMEM agree on the eta mode to ~5 significant figures; the
  # vcov agrees to ~2% (PKPDmap uses a finite-difference Jacobian over an ODE
  # solution, NONMEM uses the analytic ADVAN1 sensitivities).
  model <- PKPDsim::new_ode_model(
    code = "dAdt[0] = -(CL/V) * A[0]",
    obs = list(cmt = 1, scale = "V"),
    dose = list(cmt = 1, bioav = 1),
    parameters = list(CL = 5, V = 50)
  )
  regimen <- PKPDsim::new_regimen(
    amt = 1000, times = 0, t_inf = 1, type = "infusion"
  )
  data <- data.frame(t = c(2, 8), y = c(18, 6))
  omega <- c(0.1, 0, 0.1)
  error <- list(prop = 0.1, add = 0.1)

  fit <- function(weights) {
    suppressMessages(
      get_map_estimates(
        model = model,
        data = data,
        parameters = list(CL = 5, V = 50),
        omega = omega,
        error = error,
        regimen = regimen,
        weights = weights,
        residuals = TRUE
      )
    )
  }

  # NONMEM 7.5.1 reference output (.phi file).
  nm_unit <- list(
    eta = c(0.424754, -0.143394),
    vcov = matrix(c(0.00544459, 0.00161879,
                    0.00161879, 0.0137097), nrow = 2)
  )
  nm_weighted <- list(
    eta = c(0.440662, -0.149905),
    vcov = matrix(c(0.00240794, 0.00288528,
                    0.00288528, 0.0133440), nrow = 2)
  )

  res_unit <- fit(NULL)
  res_weighted <- fit(c(1, 3))

  # eta mode agrees tightly.
  expect_equal(as.numeric(res_unit$fit$par), nm_unit$eta, tolerance = 1e-3)
  expect_equal(as.numeric(res_weighted$fit$par), nm_weighted$eta, tolerance = 1e-3)

  # FOCE vcov agrees to within ~2% (allow 4% headroom).
  expect_equal(unname(res_unit$vcov_full), nm_unit$vcov, tolerance = 0.04)
  expect_equal(unname(res_weighted$vcov_full), nm_weighted$vcov, tolerance = 0.04)

  # Weights must actually flow into the Hessian: higher weight -> tighter vcov.
  expect_true(det(res_weighted$vcov_full) < det(res_unit$vcov_full))
})
