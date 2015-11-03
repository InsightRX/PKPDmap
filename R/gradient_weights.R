#' Gradient weighting
#' 
#' Weighting increases from 0 to 1 by a linear gradient. At `t_full_weight` the weighting is 1.
#' @param time vector of time
#' @param t_full_weight time after which weight should be 1
#' @export
gradient_weights <- function(time = NULL, 
                             t_full_weight = 0) {
                               if(is.null(t_full_weight)) {
                                 t_full_weight <- max(time)
                               }
                               gradient <- 1/t_full_weight 
                               weights <- time*gradient
                               weights[weights > 1] <- 1
                               return(weights)
                             }
