#' Join omega blocks
#' 
#' @param om1 block 1
#' @param om2 block 2
#' @param as_triangle should output block be returned as triangle (vector) or not (matrix)?
#' 
#' @export
#' 
join_blocks <- function(om1, om2, as_triangle = TRUE) {
  n1 <- lower_triangle_mat_size(om1)
  n2 <- lower_triangle_mat_size(om2)
  om1_full <- PKPDsim::triangle_to_full(om1)
  om2_full <- PKPDsim::triangle_to_full(om2)
  p1 <- cbind(om1_full, matrix(0, nrow = nrow(om1_full), ncol = ncol(om2_full)))
  p2 <- cbind(matrix(0, nrow = nrow(om2_full), ncol = ncol(om1_full)), om2_full)
  out <- rbind(p1, p2)
  if(as_triangle) {
    return(out[upper.tri(out, diag = TRUE)])
  } else {
    return(out)
  }
}
