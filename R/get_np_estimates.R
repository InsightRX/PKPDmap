#' Get non-parametric estimates given a set of support points, model, and data
#' 
#' @param parameter_grid grid (data.frame) of support points
#' @param error error model
#' @param model `PKPDsim` model
#' @param covariates `PKPDsim` covariates object
#' @param regimen `PKPDsim` regimen
#' @param data vector of observed data
#' @param weights vector of weights passed to `get_map_estimates()`
#' @param ... passed to `get_map_estimates()`
#' 
#' @export
#' 
get_np_estimates <- function(parameter_grid = NULL, 
                             error = list(prop = 0.1, add = 0.1), 
                             model = NULL,
                             covariates = NULL,
                             regimen = NULL, 
                             data = NULL,
                             weights = NULL,
                             ...) {
  all <- c()
  like <- c()
  for(i in 1:length(parameter_grid[,1])) {
    par_tmp <- list(CL = parameter_grid[i, 1], V = parameter_grid[i, 2])
    tmp <- PKPDsim::sim(
      ode = model, 
      parameters = par_tmp, 
      regimen = regimen, 
      t_obs = data$t,
      only_obs = TRUE, 
      covariates = covariates,
      checks = FALSE,
      ...)$y
    all <- rbind(all, tmp)
    like <- c(like, get_likelihood_of_data(data = data$y, ipred = tmp, error = error, weights = weights))
  }
  tmp <- data.frame(cbind(parameter_grid, like))
  tmp$CL_tmp <- tmp[,1] * tmp$like
  tmp$V_tmp  <- tmp[,2] * tmp$like 
  CL_av <- sum(tmp$CL_tmp) / sum(tmp$like)
  V_av <- sum(tmp$V_tmp) / sum(tmp$like)
  
  list(
    parameters = list(CL_av, V_av), 
    prob = tmp
  )
}
