## generic likelihood function
ll_func_generic <- function(
  data,
  # unfortunately seems no other way to do this...
  eta01, eta02, eta03, eta04, eta05, eta06, eta07, eta08, eta09, eta10,
  eta11, eta12, eta13, eta14, eta15, eta16, eta17, eta18, eta19, eta20,
  eta21, eta22, eta23, eta24,
  parameters,
  fixed = c(),
  covariates = NULL,
  covariate_names = NULL,
  regimen = regimen,
  omega_full = omega_full,
  error = error,
  model,
  t_obs,
  sig,
  weight_prior,
  ...) {
  par <- parameters
  p <- as.list(match.call())
  for(i in seq(names(par))) {
    par[[i]] <- par[[i]] * exp(p[[(paste0("eta", sprintf("%02d", i)))]])
  }
  ipred <- model(t = data[data$evid == 0,]$t,
                 parameters = par)
  y <- data$y
  res_sd <- sqrt(error$prop^2*ipred^2 + error$add^2)
  et <- mget(objects()[grep("^eta", objects())])
  et <- as.numeric(as.character(et[et != ""]))
  omega_full <- omega_full[1:length(et), 1:length(et)]
  ofv <-   c(mvtnorm::dmvnorm(et, mean=rep(0, length(et)),
                              sigma = omega_full[1:length(et), 1:length(et)],
                              log=TRUE) * weight_prior,
             stats::dnorm(y - ipred, mean = 0, sd = res_sd, log=TRUE) * weights)
  if(verbose) {
    print(ofv)
  }
  return(-sum(ofv))
}
