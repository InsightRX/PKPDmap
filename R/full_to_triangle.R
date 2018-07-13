#' Convert full matrix to lower/upper triangle (vector)
#' 
#' @param d matrix
#' 
#' @export
full_to_triangle <- function (d) {
  return(d[lower.tri(d)])
}