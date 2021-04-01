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
#' @param tdm_init_level pre-first dose TDM level
#' @param estimate_init_level estimate the pre-dose. If FALSE, will use any TDM levels before first dose as the deterministic level for compartment initiation.
#' @param init_level_weight weight in MAP fit for init level
#' @param ruv residual variability model. Required when estimate_init_level is `TRUE` to estimate error on init level.
#' @param verbose verbosity (`TRUE` or `FALSE`)
#' @export
create_iov_object <- function(cv = list(CL = 0.1),
                              omega = c(0.1),
                              bins = c(0, 24, 48, 9999),
                              parameters = list(CL = 5),
                              tdm_init_level = NULL,
                              estimate_init_level = FALSE,
                              init_level_weight = 0.5,
                              ruv = NULL,
                              fixed = NULL,
                              n = NULL,
                              verbose = TRUE) {
  if(is.null(parameters)) {
    stop("No parameters specified.")
  }
  if(is.null(omega)) {
    stop("No omega block specified.")
  }
  om_init <- NULL
  if("TDM_INIT" %in% names(parameters)) {
    if(!is.null(tdm_init_level) && tdm_init_level != 0) {
      parameters$TDM_INIT <- tdm_init_level
      if(!estimate_init_level) {
        ## Add TDM_INIT to fixed parameters, if listed in parameters and we don't want to estimate
        fixed <- unique(c(fixed, "TDM_INIT"))
      } else {
        ## otherwise we have to increase the omega_matrix with one row, to allow estimation of TDM_INIT
        if(is.null(ruv)) {
          stop("Residual variability required when estimating init level.")
        }
        init_var <- 0
        if(!is.null(ruv$exp))  init_var <- init_var + (ruv$exp)^2 # treat as proportional
        if(!is.null(ruv$prop)) init_var <- init_var + (ruv$prop)^2
        if(!is.null(ruv$add))  init_var <- init_var + (ruv$add/tdm_init_level)^2
        om_init <- init_var / init_level_weight^2
      }
    } else {
      fixed <- unique(c(fixed, "TDM_INIT"))
    }
  }
  if(is.null(cv)) {
    if(verbose) message("No IOV specified for model.")
    # make sure all kappa parameters (if present) are included in fixed vector
    iov_par <- grep("kappa_", names(parameters))
    fixed <- unique(c(fixed, names(parameters)[iov_par]))
    if(!is.null(om_init)) omega <- join_blocks(omega, om_init, as_triangle = ifelse(class(omega) == "matrix", TRUE, FALSE))
    return(list(
      parameters = parameters,
      kappa = c(),
      omega = omega,
      fixed = fixed,
      bins = bins
    ))
  }

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
  iov_par <- grepl("^kappa_", names(parameters))
  if(! any(iov_par)) {
    stop("No `kappa` parameters seem to be defined for this model.")
  }
  iov_list <- as.list(rep(0, length(kappa)))
  names(iov_list) <- kappa

  non_iov <- parameters[!iov_par]
  non_iov <- c(non_iov[! names(non_iov) %in% fixed], non_iov[fixed]) # reorder to make sure nonfixed come first
  new_par <- c(non_iov[1:n_om], iov_list, non_iov[(n_om+1):length(non_iov)])

  if(!is.null(om_init)) om_new <- join_blocks(om_new, om_init, as_triangle = ifelse(class(omega) == "matrix", TRUE, FALSE))

  return(list(
    parameters = new_par,
    kappa = kappa,
    omega = om_new,
    fixed = fixed,
    bins = bins
  ))
}
