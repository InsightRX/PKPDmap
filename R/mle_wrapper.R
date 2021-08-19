mle_wrapper <- function(minuslogl,
                        start,
                        method = "BFGS",
                        optimizer = "optim",
                        fixed = NULL,
                        data = NULL,
                        parameters = NULL,
                        control = list(reltol = 1e-5),
                        skip_hessian = FALSE) {
  
  ## Pre-processing: get eta's from global environment
  ## TODO: find better way to do this, but leave for now
  denv <- local(environment(), data)
  argnames_in_data <- names(data)[names(data) %in% names(formals(minuslogl))]
  args_in_data <- lapply(argnames_in_data, get, env = denv)
  names(args_in_data) <- argnames_in_data 

  objective_function <- function(p) {
    do.call("minuslogl", 
            args = c(p, data))
  }

  fit <- do.call("optim",
          c(
            list(
              p = start,
              fn = objective_function,
              method = method,
              hessian = FALSE,
              gr = NULL,
              control = control
            )
          ))
  
  ## Post-processing
  ## Note: may be different for different optimizers, currently supports optim only.
  if(!skip_hessian) {
    try({
      fit$hessian <- numDeriv::hessian(objective_function, fit$par)
      fit$vcov <- solve(fit$hessian, silent=TRUE)
    })
  }
  fit$coef <- fit$par
  fit$log_likelihood <- -fit$value 
  
  fit
  
}