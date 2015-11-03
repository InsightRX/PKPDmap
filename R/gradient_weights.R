#' Gradient weighting
#' 
#' Weighting increases from 0 to 1 by a linear gradient. At `t_full_weight` the weighting is 1.
#' @param time vector of time
#' @param min_weight minimum weight, default is 0.
#' @param t_full_weight time after which weight should be 1
#' @export
gradient_weights <- function(time = NULL, 
                             min_weight = 0,
                             t_full_weight = 0) {
                               if(is.null(t_full_weight)) {
                                 t_full_weight <- max(time)
                               }
                               if(min_weight < 0 | min_weight > 1) {
                                 stop("min_weight should be between 0 and 1")
                               }
                               gradient <- min_weight + (1-min_weight)/t_full_weight 
                               weights <- time*gradient
                               weights[weights > 1] <- 1
                               return(weights)
                             }
