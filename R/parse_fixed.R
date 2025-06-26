#' Parse fixed parameters
#' 
#' @inheritParams get_map_estimates
#'
parse_fixed <- function(fixed, parameters) {
  if(is.null(fixed)) return()
  fixed <- unique(fixed)
  if(length(intersect(fixed, names(parameters))) != length(fixed)) {
    warning("Warning: not all fixed parameters were found in parameter set!\n")
  }
  fixed <- names(parameters)[names(parameters) %in% fixed]
  if(length(fixed) == 0) {
    fixed <- NULL
  }
  fixed
}