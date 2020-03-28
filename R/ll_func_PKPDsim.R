#' Likelihood function for MAP optimization using PKPDsim
#' 
#' @param data vector of
#' @param sim_object design (event-table) obtained from PKPDsim to be used in simulations/optimizations
#' @param parameters parameter list
#' @param nonfixed non-fixed (i.e. estimated) parameters
#' @param error error model to use, e.g. `list(add = .5, prop = .15)`
#' @param model PKPDsim model
#' @param omega_full full omega matrix
#' @param sig signficance, used in optimization
#' @param weights weights for data, generally a vector of weights between 0 and 1.
#' @param weight_prior weight of prior
#' @param transf transformation function for data, if needed. E.g. `log(x)`. Default is `function(x) = x`.
#' @param as_eta implement as regular eta instead of exponential eta, can be vector.
#' @param censoring_idx censoring indices, used for <LOQ data.
#' @param censoring_label censoring label, used for <LOQ data
#' @param t_init init time for simulations, default 0.
#' @param iov_bins IOV bins object
#' @param calc_ofv function to calculate OFV based on simulated data and set of parameters and omega matrix
#' @param include_omega include omega in calculation of OFV?
#' @param include_error include residual error in calculation of OFV?
#' @param verbose verbose output?
#' @param eta1 eta1 
#' @param eta2 eta2 
#' @param eta3 eta3 
#' @param eta4 eta4 
#' @param eta5 eta5 
#' @param eta6 eta6 
#' @param eta7 eta7 
#' @param eta8 eta8 
#' @param eta9 eta9 
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
#' @export
ll_func_PKPDsim <- function(
  data,
  sim_object,
  parameters,
  nonfixed,
  # unfortunately seems no other way to do this using the optim package...
  eta01, eta02, eta03, eta04, eta05, eta06, eta07, eta08, eta09, eta10,
  eta11, eta12, eta13, eta14, eta15, eta16, eta17, eta18, eta19, eta20,
  eta21, eta22, eta23, eta24,
  error = error,
  model,
  omega_full,
  sig,
  weights,
  weight_prior,
  as_eta,
  censoring_idx,
  censoring_label,
  iov_bins,
  t_init = 0,
  calc_ofv,
  steady_state,
  include_omega,
  include_error,
  verbose = FALSE,
  ...) {
  par <- parameters
  p <- as.list(match.call())
  for(i in seq(nonfixed)) {
    key <- nonfixed[i]
    if(key %in% p$as_eta) {
      par[[key]] <- p[[(paste0("eta", sprintf("%02d", i)))]]
    } else {
      par[[key]] <- par[[key]] * exp(p[[(paste0("eta", sprintf("%02d", i)))]])
    }
  }
  sim_object$p <- par
  if(!is.null(steady_state)) {
    dose <- sim_object$design$dose[1]
    interval <- diff(sim_object$design[sim_object$design$evid == 1,]$t)[1]
    sim_object$A_init <- PKPDsim::calc_ss_lin(
      f = steady_state$f,
      dose,
      interval,
      model,
      parameters = par,
      covariates = covariates,
      auc = ifelse(!is.null(steady_state$auc), steady_state$auc, FALSE)
    )
    # sim_object$t_obs <- sim_object$t_obs - min(sim_object$design$t)
    # sim_object$design$t <- sim_object$design$t - min(sim_object$design$t)
  }
  ipred <- transf(PKPDsim::sim_core(
    sim_object,
    ode = model,
    duplicate_t_obs = TRUE,
    t_init = t_init)$y)
  dv <- transf(data$y)
  obs_type <- data$obs_type
  ofv_cens <- NULL
  if(!is.null(censoring_idx)) {
    cens <- data[[censoring_label]][censoring_idx] # if <0 then LLOQ, if >0 ULOQ
    cens <- ifelse(cens > 0, -1, 1)
    ipred_cens <- ipred[censoring_idx]
    ipred <- ipred[!censoring_idx]
    obs_type_cens <- obs_type[censoring_idx]
    obs_type <- obs_type[!censoring_idx]
    dv_cens <- dv[censoring_idx]
    dv <- dv[!censoring_idx]
    weights_cens <- weights[censoring_idx]
    weights <- weights[!censoring_idx]
    res_sd_cens <- sqrt(error$prop[obs_type_cens]^2*ipred_cens^2 + error$add[obs_type_cens]^2)
    ofv_cens <- stats::pnorm((dv_cens - ipred_cens) * cens, 0, res_sd_cens, log=TRUE) * weights_cens
  }
  res_sd <- sqrt(error$prop[obs_type]^2*ipred^2 + error$add[obs_type]^2)
  et <- mget(objects()[grep("^eta", objects())])
  et <- as.numeric(as.character(et[et != ""]))
  omega_full <- as.matrix(omega_full)[1:length(et), 1:length(et)]
  ofv <- calc_ofv(
    eta = et,
    omega = omega_full,
    dv = dv,
    ipred = ipred,
    res_sd = res_sd,
    weights = weights,
    weight_prior = weight_prior,
    include_omega = include_omega,
    include_error = include_error)
  ofv <- c(ofv, ofv_cens)
  if(verbose) {
    cat("-------------------------------------------------------------\n")
    cat(paste0("Eta\t: [", paste(signif(et,5), collapse=", "),"]\n"))
    cat(paste0("y_hat\t: [", paste(ipred, collapse=", "),"]\n"))
    cat(paste0("P(y)\t: [", paste(signif(exp(ofv[-1]),5), collapse=", "),"]\n"))
    cat(paste0("OFV\t: [", paste(signif(-2*sum(ofv),5), collapse=", "), "]\n"))
  }
  return(-2 * sum(ofv))
}