#' Calculate objective function value for LS fit
#'
#' @param dv observed data points
#' @param ipred individual predictions
#' @param res_sd standard deviation of the observations
#' @param weights weights vector for the observed data
#' @param ... dummy parameters, not used
#' @export
calc_ofv_ls <- function(
  dv, 
  ipred, 
  res_sd,
  weights = 1, 
  ...) {
  
  (log2pi +
    log(res_sd^2) +
    ((dv - ipred)^2 / res_sd^2)) * weights
  
}
