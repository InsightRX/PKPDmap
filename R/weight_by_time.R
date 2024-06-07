#' Gradient weighting by time
#'
#' Weighting increases from 0 to 1 by a linear gradient. At `t_end_gradient` 
#' the weighting is 1.
#' @param time vector of time
#' @param min_weight minimum weight, default is 0.
#' @param t_start_gradient time after which gradient should increase from 
#' `min_weight`. Default is 0.
#' @param t_end_gradient time after which weight should be 1. Default is 0
#'(all weights are 1)
#'
#' @export
#'
weight_by_time   <- function(time = NULL,
                             min_weight = 0,
                             t_end_gradient = NULL,
                             t_start_gradient = NULL) {
  if (is.null(t_end_gradient)) {
    t_end_gradient <- max(time)
  }
  if (is.null(t_start_gradient)) {
    t_start_gradient <- min(time)
  }
  if (min_weight < 0 |
      min_weight > 1) {
    stop("min_weight should be between 0 and 1")
  }
  if (t_end_gradient > t_start_gradient) {
    gradient <- min_weight + (1 - min_weight) / (t_end_gradient - t_start_gradient)
    weights <- rep(min_weight, length(time))
    weights[time >= t_start_gradient] <- min_weight + gradient * (time[time >= t_start_gradient] -
                                                                    t_start_gradient)
    weights[weights >= 1] <- 1
  } else {
    gradient <- 0
    weights  <- rep(1, length(time))
  }
  weights
}
