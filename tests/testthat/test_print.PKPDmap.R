test_that("MAP fit produces nice fit report", {
  
  ## Test 2cmt model but fix V2
  model2cmt <- PKPDsim::new_ode_model("pk_2cmt_iv")
  dat <- read.table(
    file = paste0(system.file(package="PKPDmap"), "/pktab1"), 
    skip=1, 
    header=TRUE
  )
  colnames(dat)[1:3] <- c("id", "t", "y")
  par2 <- list(CL = 7.67, V = 97.7, Q = 3, V2 = 50)
  tmp1 <- get_map_estimates(
    parameters = par2,
    model = model2cmt,
    regimen = PKPDsim::new_regimen(
      amt = 100000,
      times = c(0, 24),
      type = "bolus"
    ),
    omega = c(0.0406,
              0.0623, 0.117,
              0.001, 0.01, 0.1),
    fixed = c("V2"),
    error = list(prop = 0, add = sqrt(1.73E+04)),
    int_step_size = 0.1,
    data = dat[dat$id == 1, ]
  )
  
  verify_output(
    test_path("output", "print.map_estimates_1.txt"),
    PKPDmap:::print.map_estimates(tmp1)
  )
})

test_that("single value correct eta plot", {
  expect_equal(plot_eta(-2),    "o----|-----")
  expect_equal(plot_eta(-1),    "--o--|-----")
  expect_equal(plot_eta(-0.001), "----o|-----")
  expect_equal(plot_eta(0),      "-----o-----")
  expect_equal(plot_eta(+0.001), "-----|o----")
  expect_equal(plot_eta(1),      "-----|--o--")
  expect_equal(plot_eta(2),      "-----|----o")
})

test_that("vector of eta's correct plot", {
  expect_equal(
    plot_eta(c(5,1,0,-3)),
    c("-----|-----o", "-----|--o--", "-----o-----", "o-----|-----")
  )
})
