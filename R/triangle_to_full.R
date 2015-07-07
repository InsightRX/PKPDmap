#' @export
triangle_to_full <- function (vect) {
  nr <- lower_triangle_size (vect)
  k_given_i_j <- function(x , y ) ifelse( y<x, x*(x-1)/2 + y, y*(y-1)/2 + x )
  k_mat <- function(p) outer( 1:p, 1:p, k_given_i_j )
  return (matrix(vect[ k_mat( nr ) ] , nr = nr ))
}

#' @export
lower_triangle_size <- function (mat) {
  x <- length(mat)
  i <- 1
  while (x > 0) {
    x <- x-i
    i <- i+1
  }
  return(i-1)
}

#' @export
lower_triangle_vector_length <- function(comp) {
  a <- 0
  for (i in 1:comp) {
    a <- a + i
  }
  a
}