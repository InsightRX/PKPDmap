#' Calculate objective function value for MAP Bayesian fit
#'
#' @param eta eta's
#' @param omega full omega matrix
#' @param dv observed data points
#' @param ipred individual predictions
#' @param res_sd standard deviation of the observations
#' @param weights weights vector for the observed data
#' @param weight_prior relative weight of the population priors
#' @param include_omega Include deviation from population parameters in the likelihood?
#' @param include_error Include residual error in the likelihood?
#' @export
calc_ofv_map <- function(
  eta, omega,
  dv, ipred, res_sd,
  weights = 1,
  weight_prior = 1, include_omega = TRUE, include_error = TRUE) {
  ofv <-   c(mvtnorm::dmvnorm(eta, mean=rep(0, length(eta)),
                              sigma = as.matrix(omega) * 1/weight_prior,
                              log=TRUE) * include_omega,
             stats::dnorm((dv - ipred) * include_error, mean = 0, sd = res_sd, log=TRUE) * weights)
}