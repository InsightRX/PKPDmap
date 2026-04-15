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
#' @param nonfixed character vector of non-fixed parameter names
#' @param as_eta character vector of parameters estimated directly as eta
#' @param steady_state_analytic steady state settings (or NULL)
#'
calc_residuals <- function(
  obj,
  data,
  model,
  parameters_population,
  covariates,
  regimen,
  lagtime,
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
  nonfixed = NULL,
  as_eta = c(),
  steady_state_analytic = NULL,
  weight_prior_var = 1,
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
      lagtime = lagtime,
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
      lagtime = lagtime,
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

  ## Compute proper CWRES and FOCE vcov using the Jacobian (Hooker et al. 2007).
  ## The Jacobian F = df/deta is computed once via finite differences (2*n_eta
  ## simulations). Both CWRES and the FOCE Hessian/vcov are derived from F,
  ## instead of the more expensive numDeriv::hessian computation.
  if (!is.null(nonfixed) && length(nonfixed) > 0) {
    n_before <- nrow(data_before_init)
    # Extract ipred for real observations only (excluding before-init)
    ipred_real <- sim_ipred$y[(n_before + 1):length(sim_ipred$y)]

    foce_result <- tryCatch(
      calc_cwres(
        eta_hat = obj$fit$coef,
        ipred_raw = ipred_real,
        y = data$y,
        omega_full = omega_full,
        error = error,
        obs_type = data$obs_type,
        transf = transf,
        model = model,
        parameters_population = parameters_population,
        nonfixed = nonfixed,
        as_eta = as_eta,
        covariates = covariates,
        regimen = regimen,
        lagtime = lagtime,
        t_obs = data$t,
        obs_type_sim = data$obs_type,
        int_step_size = int_step_size,
        iov_bins = iov_bins,
        A_init = A_init_individual,
        t_init = t_init,
        steady_state_analytic = steady_state_analytic,
        weight_prior_var = weight_prior_var,
        ...
      ),
      error = function(e) {
        warning("CWRES computation failed: ", e$message)
        list(cwres = rep(NA_real_, length(data$y)), vcov = NULL)
      }
    )

    # Pad with zeros for before-init observations and apply weights
    obj$cwres <- c(rep(0, n_before), foce_result$cwres * weights)

    # Set censored observations to NA
    if (!is.null(censoring) && any(censoring_idx)) {
      obj$cwres[censoring_idx] <- NA_real_
    }

    # Store FOCE vcov for use in get_map_estimates
    obj$foce_vcov <- foce_result$vcov
  }

  ## Add covariates and parameters to obj
  if(output_include$covariates && !is.null(covariates)) {
    obj$covariates_time <- sim_ipred[!duplicated(sim_ipred$t), names(covariates)]
  }
  if(output_include$parameters) {
    obj$parameters_time <- sim_ipred[!duplicated(sim_ipred$t), names(parameters)]
  }

  obj
}