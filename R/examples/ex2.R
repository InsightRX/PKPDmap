system.time({
  library(bbmle)
  library(mvtnorm)
  library(PKPDsim)
  library(PKPDmap)
  library(pk1cmtiv)
  model <- pk1cmtiv::model("pk_1cmt_iv")
  regimen <- new_regimen(amt = 1000, interval=12, t_inf= 1, n = 6, type="infusion")
  data <- data.frame(cbind(t=23.75, y=12, evid=0))
  omega <- cv_to_omega(list(CL = 0.2816, Vc = 0.3715), parameters = list(CL = 1.08, Vc = 68))
})

pk_iv_1cmt_R <- function (t, t_inf = 1, tau = 24, dose=120, 
                          parameters = list (CL = 0.345, V = 1.75), 
                          ruv = NULL) {
  # takes either t as a vector or CL/Vc, but not both
  tmp <- with(parameters, {
    k <- CL / V
    tmp <- c()
    tmp <- c(tmp, (dose / (CL * t_inf)) * (1-exp(-k*t[t < t_inf])) )
    tmp <- c(tmp, (dose / (CL * t_inf)) * (1-exp(-k*t_inf)) * exp(-k*(t[t >= t_inf] - t_inf)) )
    if(!is.null(ruv)) {
      tmp <- add_ruv (tmp, ruv)
    }
    tmp
  })
  tmp
}

plot(log(pk_iv_1cmt_R(1:100)))

system.time({
  fit <- get_map_estimates (
    model = pk1cmtiv::model("pk_1cmt_iv"), 
    int_step_size = .1,
#    model = pk_1cmt_iv_R,
    data = data, 
    parameters = list(CL=5.289062, V=68.6),
    omega = omega,
    regimen = regimen,
    error = list(prop = 0, add = 3.5),
    method= "BFGS"
  )
})
