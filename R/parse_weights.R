#' Parse weights argument
#' 
#' @inheritParams get_map_estimates
#' 
parse_weights <- function(weights, t_obs) {
  if(!is.null(weights)) {
    if(length(weights) != length(t_obs)) {
      stop("Vector of weights of different size than observation vector!")
    }
  } else {
    weights <- rep(1, length(t_obs))
  }
  weights
}