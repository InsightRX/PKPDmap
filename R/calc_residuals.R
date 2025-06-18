#' Function to calculate residuals from a fit,
#' based on a fitted parameter set, and data object.
#' 
#' @inheritParams get_map_estimates
#' 
calc_residuals <- function(
  obj,
  data,
  model,
  parameters_population,
  covariates,
  regimen,
  omega_full,
  error,
  weights,
  transf,
  t_obs,
  obs_type,
  A_init_population,
  A_init_individual,
  t_init = 0,
  iov_bins = NULL,
  output_include = c(),
  int_step_size = 0.01,
  censoring = NULL,
  censoring_idx = NULL,
  data_before_init,
  ltbs = FALSE,
  ...
) {
  suppressMessages({
    ## After fitting individual parameters, don't pass the mixture group to the simulator 
    ## (so `mixture=NULL`), otherwise `sim()` will use the population value for the 
    ## specified group, and not the individual fitted parameter.
    sim_ipred <- PKPDsim::sim_ode(
      ode = model,
      parameters = obj$parameters,
      mixture_group = NULL,
      covariates = covariates,
      n_ind = 1,
      int_step_size = int_step_size,
      regimen = regimen,
      t_obs = t_obs,
      obs_type = obs_type,
      only_obs = TRUE,
      checks = FALSE,
      A_init = A_init_individual,
      iov_bins = iov_bins,
      output_include = output_include,
      t_init = t_init,
      ...
    )
  })
  suppressMessages({
    sim_pred <- PKPDsim::sim_ode(
      ode = model,
      parameters = parameters_population,
      mixture_group = NULL,
      covariates = covariates,
      n_ind = 1,
      int_step_size = int_step_size,
      regimen = regimen,
      t_obs = c(data_before_init$t, t_obs),
      obs_type = c(data_before_init$obs_type, data$obs_type),
      only_obs = TRUE,
      checks = FALSE,
      iov_bins = iov_bins,
      A_init = A_init_population,
      t_init = t_init,
      ...
    )
  })
  ipred <- sim_ipred$y
  pred <- sim_pred$y
  w_ipred <- sqrt(error$prop[data$obs_type]^2 * transf(ipred)^2 + error$add[data$obs_type]^2)
  w_pred <- sqrt(error$prop[data$obs_type]^2 * transf(pred)^2 + error$add[data$obs_type]^2)
  if(!(all(data$t == sim_ipred$t) && all(data$obs_type == sim_ipred$obs_type))) {
    warning("Mismatch in times and observation typese between input data and predictions. Be careful interpreting results from fit.")
  }
  y <- data$y
  cf <- obj$fit$coef
  prob <- list(par = c(mvtnorm::pmvnorm(cf, mean=rep(0, length(cf)),
                                        sigma = omega_full[1:length(cf), 1:length(cf)])),
               data = stats::pnorm(transf(y) - transf(ipred), mean = 0, sd = w_ipred))
  obj$res <- (transf(y) - transf(pred))
  obj$weights <- c(rep(0, length(data_before_init$t)), weights)
  obj$wres <- (obj$res / w_pred) * obj$weights
  obj$cwres <- obj$res / sqrt(abs(cov(transf(pred), transf(y)))) * c(rep(0, nrow(data_before_init)), obj$weights)
  # Note: in NONMEM CWRES is on the population level, so can't really compare. NONMEM calls this CIWRES, it seems.
  obj$ires <- (transf(y) - transf(ipred))
  obj$iwres <- (obj$ires / w_ipred)
  if(is.null(censoring)) {
    obj$censoring <- rep(0, length(y))
  } else {
    obj$censoring <- data[[censoring]]
    if(any(censoring_idx)) { # turn probabilities into "IWRES"-equivalent for censored data
      obj$iwres[censoring_idx] <- calc_res_from_prob(prob$data[censoring_idx])
      obj$ires[censoring_idx] <- obj$iwres[censoring_idx] * w_ipred[censoring_idx]
      ## if we would calculate the likelihood for the data population parameters given the data, 
      ## we could also calculate the the equivalents for cwres, wres, and res. However we 
      ## currently don't have a need to calculate. So setting to NA to avoid wrong interpretation.
      obj$cwres[censoring_idx] <- NA_real_
      obj$wres[censoring_idx] <- NA_real_
      obj$res[censoring_idx] <- NA_real_
    }
  }
  obj$iwres_weighted <- obj$iwres * obj$weights
  obj$pred <- pred
  obj$ipred <- ipred
  obj$prob <- prob
  obj$dv <- y
  obj$obs_type <- sim_ipred$obs_type
  if(output_include$covariates && !is.null(covariates)) {
    obj$covariates_time <- sim_ipred[!duplicated(sim_ipred$t), names(covariates)]
  }
  if(output_include$parameters) {
    obj$parameters_time <- sim_ipred[!duplicated(sim_ipred$t), names(parameters)]
  }
  obj$mahalanobis <- get_mahalanobis(
    y, 
    ipred, 
    w_ipred, 
    ltbs
  )
  obj
}