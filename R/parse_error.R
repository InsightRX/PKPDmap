#' Parse specified error magnitudes
#' 
#' @inheritParams get_map_estimates
#'
parse_error <- function(error) {
  if(!is.null(error)) { ## ensure all error types are specified and not NULL
    if(is.null(error$prop)) error$prop <- 0
    if(is.null(error$add)) error$add <- 0
    if(is.null(error$exp)) error$exp <- 0
    if(sum(unlist(error)) == 0) {
      stop("No residual error model specified, or residual error is 0.")
    }
  } else {
    stop("No residual error specified, cannot perform fit.")
  }
  error
}
