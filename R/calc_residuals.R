#' Function to calculate residuals from a fit,
#' based on a fitted parameter set, and data object.
#' 
#' @inheritParams get_map_estimates
#' @param obj temporary storage object in `get_map_estimates()`
#' @param parameters_population parameters for population predictions
#' @param omega_full full omega matrix
#' @param transf transformation function
#' @param A_init_population init vector for model state (population)
#' @param A_init_individual init vector for model state (individual)
#' @param censoring_idx vector with indices for censoring
#' @param data_before_init data.frame with data before initial dose
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
  A_init_population,
  A_init_individual,
  t_init = 0,
  iov_bins = NULL,
  output_include = c(),
  int_step_size = 0.01,
  censoring = NULL,
  censoring_idx = NULL,
  data_before_init = NULL,
  ltbs = FALSE,
  ...
) {

  ## Observation vectors
  t_obs <- c(data_before_init$t, data$t)
  obs_type <- c(data_before_init$obs_type, data$obs_type)
  
  ## Perform simulations for ipred and pred
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

  ## Parse the residuals from the predictions and the data, and add to obj
  obj <- parse_residuals_from_predictions(
    obj,
    sim_ipred,
    sim_pred,
    data,
    omega_full,
    transf,
    error,
    censoring,
    censoring_idx,
    data_before_init,
    weights
  )
  
  ## Add covariates and parameters to obj
  if(output_include$covariates && !is.null(covariates)) {
    obj$covariates_time <- sim_ipred[!duplicated(sim_ipred$t), names(covariates)]
  }
  if(output_include$parameters) {
    obj$parameters_time <- sim_ipred[!duplicated(sim_ipred$t), names(parameters)]
  }

  obj
}