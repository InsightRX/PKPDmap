#' Generic likelihood function
#' 
#' @inheritParams ll_func_PKPDsim
#'
#' @param omega_full full omega matrix
#' @param sig signficance, used in optimization
#' @param eta01 eta1 
#' @param eta02 eta2 
#' @param eta03 eta3 
#' @param eta04 eta4 
#' @param eta05 eta5 
#' @param eta06 eta6 
#' @param eta07 eta7 
#' @param eta08 eta8 
#' @param eta09 eta9 
#' @param eta10 eta10 
#' @param eta11 eta11 
#' @param eta12 eta12 
#' @param eta13 eta13 
#' @param eta14 eta14 
#' @param eta15 eta15 
#' @param eta16 eta16 
#' @param eta17 eta17 
#' @param eta18 eta19
#' @param eta19 eta19
#' @param eta20 eta20
#' @param eta21 eta21
#' @param eta22 eta22
#' @param eta23 eta23
#' @param eta24 eta24 
#' @param ... passed on to PKPDsim
#' 
ll_func_generic <- function(
  data,
  # unfortunately seems no other way to do this...
  eta01, eta02, eta03, eta04, eta05, eta06, eta07, eta08, eta09, eta10,
  eta11, eta12, eta13, eta14, eta15, eta16, eta17, eta18, eta19, eta20,
  eta21, eta22, eta23, eta24,
  parameters,
  covariates = NULL,
  regimen = regimen,
  lagtime = NULL,
  omega_full = omega_full,
  error = error,
  model,
  sig,
  verbose = FALSE,
  ...
) {
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
                              log=TRUE),
             stats::dnorm(y - ipred, mean = 0, sd = res_sd, log=TRUE) * weights)
  if(verbose) {
    print(ofv)
  }

  -sum(ofv)
}
