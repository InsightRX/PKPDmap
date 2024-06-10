#' Convert full matrix to lower/upper triangle (vector)
#' 
#' @param d matrix
#' 
#' @export
#' 
full_to_triangle <- function (d) {
  d[lower.tri(d)]
}