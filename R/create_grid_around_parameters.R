#' Create an expanded multivariate grid of parameters around a given set of parameters
#' 
#' The expanded grid outputted from this function can be used for performing 
#' extended grid non-parametric estimation
#' 
#' @param parameters list of parameters
#' @param grid_size number of points per parameters
#' @param span width of the grid in either side
#' @export
create_grid_around_parameters <- function(parameters = list(),
                                          grid_size = 4,
                                          span = 0.2) {
  lst <- list()
  for(i in names(parameters)) {
    lst[[i]] <- unlist(parameters)[i] * (1 + seq(from = -span, to = span, length.out = n))
  }
  return(do.call("expand.grid", lst))
}
