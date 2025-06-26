#' Wrapper for non-parametric fit
#' 
#' @param obj fit object
#' @inheritParams get_map_estimates
#' 
np_fit_wrapper <- function(
  obj,
  model,
  regimen,
  data,
  covariates = NULL,
  weights = 1,
  np_settings = list()
) {
  par <- obj$parameters
  obj$parameters_map <- par ## keep copy of MAP estimates
  np_settings <- replace_list_elements(
    np_settings_default, 
    np_settings
  )
  if(is.null(np_settings$error)) { # if no specific error magnitude is specified for NP, just use same as used for MAP
    np_settings$error <- error
  }
  if(np_settings$grid_adaptive) { # do a first pass with a broad grid
    pars_grid <- create_grid_around_parameters(
      parameters = par,
      span = np_settings$grid_span_adaptive,
      exponential = np_settings$grid_exponential_adaptive,
      grid_size = np_settings$grid_size_adaptive
    )
    np <- get_np_estimates(
      parameter_grid = pars_grid,
      error = np_settings$error,
      model = model,
      regimen = regimen,
      data = data$y,
      t_obs = data$t,
      covariates = covariates,
      weights = weights
    )
    # take the estimates with highest probability as starting point for next grid
    tmp <- np$prob[order(-np$prob$like),][1,]
    par <- as.list(tmp[1:length(np$parameters)])
  }
  pars_grid <- create_grid_around_parameters(
    parameters = par,
    span = np_settings$grid_span,
    exponential = np_settings$grid_exponential,
    grid_size = np_settings$grid_size
  )
  np <- get_np_estimates(
    parameter_grid = pars_grid,
    error = np_settings$error,
    model = model,
    regimen = regimen,
    data = data$y,
    t_obs = data$t,
    covariates = covariates,
    weights = weights
  )
  for(i in 1:length(par)) {
    par[[i]] <- np$parameters[[i]]
  }
  obj$np <- list(prob = np$prob)
  obj$parameters <- par
  obj
}