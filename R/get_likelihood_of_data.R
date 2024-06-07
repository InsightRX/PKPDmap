#' Get the likelihood of observed data given a specific error function
#' 
#' @param data vector of observed data points
#' @param ipred vector of predicted values,
#' @param error list specifying error model (list of `prop` and `add`), assuming normal distribution
#' @param weights vector of weights for observations (not used by default, `NULL`)
#' 
#' @export
#' 
get_likelihood_of_data <- function(data, ipred, error, weights = NULL) {
  if(is.null(weights)) {
    weights <- rep(1, length(ipred))
  }
  res_sd <- sqrt(error$prop^2 * ipred^2 + error$add^2)
  if(!is.null(weights) && length(weights) == length(ipred)) {
    likelhd <- stats::dnorm((data - ipred) * weights, mean = 0, sd = res_sd, log = FALSE)
  } else {
    likelhd <- stats::dnorm((data - ipred), mean = 0, sd = res_sd, log = FALSE)
  }
  prod(likelhd)
}
