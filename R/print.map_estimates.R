#' Print MAP estimates
#' 
#' @param x output object from get_map_estimates
#' @param ... etc
#' @export
print.map_estimates <- function(x, ...) {
  for (i in seq(x$parameters)) {
    cat(paste0(names(x$parameters[i]), ": ", round(as.numeric(x$parameters[i]), 3), "\n"))
  }
}