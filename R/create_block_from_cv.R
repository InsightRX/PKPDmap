#' Create an variance-covariance matrix from specified CV
#' Mostly useful for creating IOV blocks
#' 
#' @param cv CV of 
#' @param n number of occassions to allow
#' 
#' @export
#' 
create_block_from_cv <- function(cv = 0.1, n = 4) {
  diag(n)[upper.tri(diag(n), diag = TRUE)] * rep(cv^2)
}
