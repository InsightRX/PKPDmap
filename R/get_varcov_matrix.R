#' Get the variance-covariance matrix from the fit object.
#' Includes some safety checks.
#' 
#' @param obj fit object from get_map_estimates
#' @param fallback a fallback variance-covariance matrix, if the varcov-matrix
#' from the fit object does not satisfy safety checks (e.g. not postive-
#' definite).
#' 
get_varcov_matrix <- function(
    obj, 
    fallback = NULL
) {
  if(is.null(obj$vcov_full)) {
    if(!is.null(fallback)) {
      obj$vcov_full <- fallback
    } else {
      return(NULL)
    }
  } 
  if(any(is.na(obj$vcov_full)) || !PKPDsim::is_positive_definite(obj$vcov_full)) {
    if(!is.null(fallback)) {
      obj$vcov_full <- fallback
      warning("Var-cov matrix of MAP estimate not positive-definite, returning original `omega` instead.")
    }
  }
  return(obj$vcov_full[t(!upper.tri(obj$vcov_full))])
}
