#' Wrapper function around optimization function (e.g. `optim`)
#'
#' Originally the function used for this was bbmle::mle(), but that
#' function was needlessly complex and had several outstanding
#' development issues.
#'
#' @param minuslogl function that calculates -log(L) based on parameters and other additional `data` passed to it.
#' @param start initial estimates for the parameters
#' @param method estimation method passed to optimization library.
#' @param optimizer optimizer used. Currently only `optim` is supported.
#' @param data additional data required for `minuslogl` function to perform optimization.
#' @param control additional parameters that control the optimization process.
#' @param skip_hessian skip the calculation of the Hessian and variance-covariance matrix at the MAP estimates?
#' 
mle_wrapper <- function(minuslogl,
                        start,
                        method = "BFGS",
                        optimizer = c("optim"),
                        data = NULL,
                        control = list(reltol = 1e-5),
                        skip_hessian = FALSE) {

  optimizer <- match.arg(optimizer)
  
  ## get eta's from global environment
  ## TODO: find better way to do this, but leave for now
  denv <- local(environment(), data)
  argnames_in_data <- intersect(names(data), names(formals(minuslogl)))
  args_in_data <- lapply(argnames_in_data, get, env = denv)
  names(args_in_data) <- argnames_in_data
  
  objective_function <- function(par) {
    do.call("minuslogl", args = c(par, data))
  }
  
  ## Perform optimization
  fit <- do.call(
    optimizer,
    c(
      list(
        par = start,
        fn = objective_function,
        method = method,
        hessian = FALSE,
        gr = NULL,
        control = control
      )
    ))
  
  ## Post-processing
  ## Note: may be different for different optimizers, currently supports optim only.
  if (!skip_hessian) {
    ## optimization libraries are commonly also able to calculate the Hessian, 
    ## but re-implementing here wrapped in try() to make it safer.
    try(
      {
        fit$hessian <- numDeriv::hessian(objective_function, fit$par)
        fit$vcov <- solve(fit$hessian)
      }, 
      silent = FALSE
    )
  }
  fit$coef <- fit$par
  fit$log_likelihood <- -fit$value
  
  fit
  
}