#' Calculates the Mahalanobis distance from observations, predictions and residual error-weighted predictions
#' 
#' @param y observed values (DV)
#' @param ipred individualized predicted values
#' @param w_ipred weighted individualized predicted values
#' @param ltbs logical indicating whether both sides should be log transformed
#'
get_mahalanobis <- function(y, ipred, w_ipred, ltbs = FALSE){
  
  if (length(w_ipred) == 1) {
    use_cov <- w_ipred ** 2
  } else {
    use_cov <- diag(w_ipred ** 2)
  }
  
  if (ltbs) {
    y <- log(y)
    ipred <- log(ipred)
  }
  
  tryCatch(
    expr = {
      stats::mahalanobis(y, ipred, cov = use_cov)
    },
    error = function(e){ 
      NULL
    }
  )
}