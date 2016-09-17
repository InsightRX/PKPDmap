#' Get the likelihood of observed data given a specific error function
#' 
#' @param data vector of observed data points
#' @param ipred vector of predicted values,
#' @param error list specifying error model (list of `prop` and `add`), assuming normal distribution
#' @export
get_likelihood_of_data <- function(data, ipred, error) {
  res_sd <- sqrt(error$prop^2 * ipred^2 + error$add^2)
  likelhd <- prod(stats::dnorm((data - ipred), mean = 0, sd = res_sd, log = FALSE))
  return(likelhd)
}
