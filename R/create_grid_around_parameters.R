#' Create an expanded multivariate grid of parameters around a given set of parameters
#' 
#' The expanded grid outputted from this function can be used for performing 
#' extended grid non-parametric estimation
#' 
#' @param parameters list of parameters
#' @param grid_size number of points per parameters
#' @param exponential use exponential grid?
#' @param span width of the grid in either side
#' @param fix vector of fixed parameter names
#' 
#' @export
#' 
create_grid_around_parameters <- function(parameters = list(),
                                          grid_size = 4,
                                          exponential = FALSE,
                                          span = 0.5,
                                          fix = NULL) {
  lst <- list()
  if(!is.null(fix)) {
    for(i in fix) {
      parameters[[fix]] <- NULL
    }
  }
  if(length(parameters) < 2) {
    stop("At least 2 unfixed parameters are expected to be estimated.")
  }
  if(exponential) {
    for(i in names(parameters)) {
      lst[[i]] <- unlist(parameters)[i] * exp(seq(from = -span, to = span, length.out = grid_size))
    }
  } else {
    if(span >= 1) {
      stop("Span cannot be >1 with non-exponential grid.")
    }
    for(i in names(parameters)) {
      lst[[i]] <- unlist(parameters)[i] * (1 + seq(from = -span, to = span, length.out = grid_size))
    }
  }
  
  do.call("expand.grid", lst)
}
