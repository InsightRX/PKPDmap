#' Parse the weight of the prior, make sure
#' it is on the variance scale.
#' @param weight_prior weight of the prior specified on SD scale
#' @param type type of fit, e.g. `"map"` or `"pls"`
#'  
parse_weight_prior <- function(
  weight_prior = 1,
  type = "map"
) {
  if(is.null(weight_prior) || is.na(weight_prior)) {
    weight_prior <- 1
  }
  weight_prior_var <- weight_prior^2
  if(tolower(type) == "pls") {
    ## RK: Empirically determined to be a good penalty.
    ##     In principle, PLS is just MAP with very flat priors
    weight_prior_var <- 0.001
  }
  weight_prior_var
}