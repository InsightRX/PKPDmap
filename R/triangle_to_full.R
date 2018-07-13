#' Convert lower/upper triangle of a matrix to full matrix
#' 
#' @param vect vector describing the lower or upper triangle of a matrix
#' 
#' @export
triangle_to_full <- function (vect) {
  nr <- lower_triangle_mat_size (vect)
  k_given_i_j <- function(x , y ) ifelse( y<x, x*(x-1)/2 + y, y*(y-1)/2 + x )
  k_mat <- function(p) outer( 1:p, 1:p, k_given_i_j )
  return (matrix(vect[ k_mat( nr ) ] , nrow = nr ))
}