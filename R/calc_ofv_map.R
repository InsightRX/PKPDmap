#' Calculate objective function value for MAP Bayesian fit
#'
#' @param eta eta's specified as row matrix, e.g. `matrix(c(...), nrow=1)`.
#' @param omega full omega matrix
#' @param omega_inv inverse of omega matrix. Passed to this function to avoid doing computations in each iteration of the search.
#' @param omega_eigen eigenvalue decomposation (logged) of omega matrix. Passed to this function to avoid doing computations in each iteration of the search.
#' @param dv observed data points
#' @param ipred individual predictions
#' @param res_sd standard deviation of the observations
#' @param weights weights vector for the observed data
#' @param include_omega Include deviation from population parameters in the likelihood?
#' @param include_error Include residual error in the likelihood?
#' @export
calc_ofv_map <- function(
  eta, 
  omega,
  omega_inv,
  omega_eigen,
  dv, 
  ipred, 
  res_sd,
  weights = 1,
  include_omega = TRUE, 
  include_error = TRUE) {
  
  # The code below implements:
  #
  # ofv <-  -2 * c(
  #   include_omega * mvtnorm::dmvnorm(
  #     eta, 
  #     mean = rep(0, length(eta)),
  #     sigma = as.matrix(omega),
  #     log = TRUE
  #   ),
  #   weights * stats::dnorm(
  #     (dv - ipred) * include_error, 
  #     mean = 0, 
  #     sd = res_sd, 
  #     log=TRUE
  #   )
  # )
  
  ofv_om <- log2pi * ncol(omega) +
            omega_eigen +
            include_omega * diag(eta %*% omega_inv %*% t(eta))
  ofv_y <- (log2pi +
            log(res_sd^2) +
            ((dv - ipred)^2 * include_error / res_sd^2) * weights)
  
  c(ofv_om, ofv_y)
  
}
