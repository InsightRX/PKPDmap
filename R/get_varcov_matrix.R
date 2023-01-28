#' Get the variance-covariance matrix from the fit object.
#' Includes some safety checks.
#' 
#' @param obj fit object from get_map_estimates
#' @param fallback a fallback variance-covariance matrix, if the varcov-matrix
#' from the fit object does not satisfy safety checks (e.g. not postive-
#' definite).
#' 
get_varcov_matrix <- function(
    vcov, 
    fallback = NULL
) {
  vcov_full <- vcov
  if(is.null(vcov_full)) {
    if(!is.null(fallback)) {
      vcov_full <- fallback
    } else {
      return(NULL)
    }
  } 
  if(any(is.na(vcov_full)) || !PKPDsim::is_positive_definite(vcov_full)) {
    if(!is.null(fallback)) {
      vcov_full <- fallback
      warning("Var-cov matrix of MAP estimate not positive-definite, returning original `omega` instead.")
    } else {
      warning("Var-cov matrix of MAP estimate not positive-definite, no fallback specified.")
      return(NULL)
    }
  }
  return(vcov_full)
}
