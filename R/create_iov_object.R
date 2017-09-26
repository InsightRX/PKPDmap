#' Calculate the specs for including IOV, e.g. which elements should be fixed
#' and the right size of the new omega matrix.
#' Currently only working for a single parameter!
#'
#' @param cv list of CVs of parameters
#' @param omega omega matrix (lower triangle)
#' @param bins vector of bins for IOV
#' @param parameters named list of parameter values
#' @param fixed vector of fixed parameters
#' @param n number of IOV elements, will be determined automatically from data, `n` will override if not NULL.
#' @param reverse_time reverse the timing of the bins to mean "t before last observation". `FALSE` by default.
#' @param verbose verbosity (`TRUE` or `FALSE`)
#' @export
create_iov_object <- function(cv = list(CL = 0.1),
                              omega = c(0.1),
                              bins = c(0, 24, 48, 9999),
                              parameters = list(CL = 5),
                              fixed = NULL,
                              n = NULL,
                              verbose = TRUE) {
  if(is.null(parameters)) {
    stop("No parameters specified.")
  }
  if(is.null(omega)) {
    stop("No omega block specified.")
  }
  if(is.null(cv)) {
    if(verbose) message("No IOV specified for model.")
    # make sure all kappa parameters (if present) are included in fixed vector
    iov_par <- grep("kappa_", names(parameters))
    fixed <- unique(c(fixed, names(parameters)[iov_par]))
    return(list(
      parameters = parameters,
      kappa = c(),
      omega = omega,
      fixed = fixed,
      bins = bins
    ))
  }

  # for(i in seq(parameters)) {
  #
  # }
  ## 1. get n from data / bins
  if(is.null(n)) n <- length(bins) - 1
  kappa <- c()
  om_new <- omega
  for(i in seq(names(cv))) {
    kappa <- c(kappa, paste0("kappa_", names(cv)[i], "_", 1:n))
    om2 <- create_block_from_cv(cv = cv[[i]], n)
    om_new <- join_blocks(om_new, om2)
  }
  n_om <- lower_triangle_mat_size(omega)

  ## reshuffle parameters
  iov_par <- !is.na(c(stringr::str_match(names(parameters), "kappa_")))
  if(! any(iov_par)) {
    stop("No `kappa` parameters seem to be defined for this model.")
  }
  iov_list <- as.list(rep(0, length(kappa)))
  names(iov_list) <- kappa

  non_iov <- parameters[!iov_par]
  new_par <- c(non_iov[1:n_om], iov_list, non_iov[(n_om+1):length(non_iov)])
  return(list(
    parameters = new_par,
    kappa = kappa,
    omega = om_new,
    fixed = fixed,
    bins = bins
  ))
}
