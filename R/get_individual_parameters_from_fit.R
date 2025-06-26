#' Get the individual parameters from the fitted coefficients (etas)
#' 
#' @param fit fit object
#' @param parameters list of population parameters
#'
get_individual_parameters_from_fit <- function(
  fit,
  parameters,
  nonfixed = NULL,
  as_eta = NULL
) {
  par <- parameters
  for(i in seq(nonfixed)) {
    key <- nonfixed[i]
    if(key %in% as_eta) {
      par[[key]] <- as.numeric(fit$coef[i])
    } else {
      par[[key]] <- as.numeric(as.numeric(par[[key]]) * exp(as.numeric(fit$coef[i])))
    }
  }
  par
}
