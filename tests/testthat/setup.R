mod_1cmt_oral_lagtime <- PKPDsim::new_ode_model(
  code = "
    dAdt[0] = -KA * A[0]
    dAdt[1] = +KA * A[0] -(CL/V) * A[1]
  ",
  lagtime = c("TLAG", 0),
  obs = list(cmt = 2, scale = "V"),
  dose = list(cmt = 1, bioav = 1),
  parameters = list(CL = 5, V = 50, KA = 0.5, TLAG = 0.83)
)