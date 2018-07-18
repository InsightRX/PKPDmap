#' Get the size of matrix specified as vector elements of the lower triangle
#' 
#' @param mat vector of matrix elements of lower triangle
#' @export
lower_triangle_mat_size <- function (mat) {
  x <- length(mat)
  i <- 1
  while (x > 0) {
    x <- x-i
    i <- i+1
  }
  return(i-1)
}
