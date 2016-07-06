#' Iterative Two State (ITS) estimation method for PK data
#' 
#' Note: this estimation method is mostly for demonstration purposes
#' as it is unstable, biased, and imprecise. The implementation
#' in this package is also pretty slow. It is however
#' a great example estimation method e.g. for teaching the concept
#' of estimation methods in pharmacokinetics.
#'       
#' @param parameters list of intial parameter estimates
#' @param omega omega lower triangle (IIV)
#' @param err additive error
#' @param regimen `PKPDsim` regimen
#' @param model `PKPDsim` model
#' @param max_iter maximum number of iterations
#' @param data observed data
#' @param ... additional arguments passed to `get_map_estimates` function
#' @export
run_its <- function(
  parameters = NULL, 
  omega = NULL, 
  err = 100,
  regimen = NULL,
  model = NULL, 
  max_iter = 5, 
  data = NULL,
  ...) {
  par_tmp <- parameters
  omega_tmp <- omega
  err <- 100
  for(j in 1:max_iter) {
    pars <- c()
    etas <- c()
    message(paste0("ITS iteration: ", j))
    for(i in seq(unique(dat$id))) {
      dv <- data[data$id == i & dat$EVID == 0,]
      tmp <- get_map_estimates(parameters = par_tmp,
                               model = model,
                               regimen = regimen,
                               omega = omega_tmp,
                               weights = rep(1, length(dv[,1])),
                               error = list("add" = err, prop=0),
                               data = dv,
                               residuals = TRUE,
                               ...
      )
      pars <- rbind(pars, cbind(t(as.numeric(unlist(tmp$parameters)))))
      etas <- rbind(etas, cbind(t(as.numeric(unlist(coef(tmp$fit))))))
    }
    par_mn <- apply(pars, 2, mean)
    par_tmp <- parameters
    for(k in seq(names(parameters))) {
      par_tmp[[k]] <- par_mn[k]
    }
    omega_tmp <- cov(etas)
    cvs <- diag(omega_tmp)
    omega_tmp <- omega_tmp[!lower.tri(omega_tmp)]
    err <- sd(tmp$res)
    message(" - Population means: ", paste(round(par_mn, 3), collapse = "  "))
    message(" - IIV %: ", paste(round(100 * sqrt(cvs), 1), collapse = "  "))
    message(" - Error: ", err)
    if(any(round(100 * sqrt(cvs), 1) == 0)) {
      stop("IIV zero for some parameter, stopping search.")
    }
  }
} 
