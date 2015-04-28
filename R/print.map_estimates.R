#' @export
print.map_estimates <- function(obj) {
  for (i in seq(obj$parameters)) {
    cat(paste0(names(obj$parameters[i]), ": ", round(as.numeric(obj$parameters[i]), 3), "\n"))
  }
}