#' Calculate objective function value for LS fit
#'
#' The OFV is defined as -2 times the log of the density of the normal distribution 
#' at the point `dv-ipred` with mean 0 and standard deviation res_sd.
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
