#' Parse omega matrix into full matrix, and
#' perform some checks
#' 
#' @inheritParams get_map_estimates
#' 
parse_omega_matrix <- function(
  omega,
  parameters,
  fixed
) {
  if(inherits(omega, "matrix")) {
    omega_full <- omega # dummy om matrix
  } else {
    omega_full <-  PKPDsim::triangle_to_full(omega) # dummy om matrix
  }
  om_nonfixed <-  PKPDsim::triangle_to_full(omega)
  if(nrow(om_nonfixed) < (length(parameters) - length(fixed))) {
    msg <- "Provided omega matrix is smaller than expected based on the number of model parameters. Either fix some parameters or increase the size of the omega matrix.\n"
    msg <- c(msg,
             paste0("Non-fixed omegas: ", paste(om_nonfixed, collapse=", "), "\n"),
             paste0("Parameters: ", paste(parameters, collapse=", "), "\n"),
             paste0("Fixed: ", paste(fixed, collapse=","), "\n"))
    stop(msg)
  }
  omega_full
}