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
  nonfixed <- names(parameters)[is.na(match(names(parameters), fixed))]
  n_nonfix <- length(nonfixed)
  eta <- list()
  for(i in seq(nonfixed)) {
    eta[[paste0("eta", sprintf("%02d", i))]] <- 0
  }
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
  omega_full_est <- omega_full[
    1:n_nonfix, 
    1:n_nonfix
  ]
  list(
    full = omega_full,    # full omega block
    est = omega_full_est, # full omega block but only for non-fixed parameters,
    fixed = fixed,        # vector of fixed parameters
    nonfixed = nonfixed,  # vector of non-fixed (estimaterd) parameters
    eta = eta             # list of etas
  )
}