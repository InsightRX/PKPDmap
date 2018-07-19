#' Crude implementation of Iterative Two State Bayesian (ITSB) estimation method for PopPK data
#' 
#' Note: this estimation method is mostly for demonstration purposes
#' as it is unstable, biased, and imprecise. The implementation
#' in this package is also pretty slow. It is however
#' a nice example estimation method e.g. for teaching the concept
#' of estimation methods in pharmacokinetics.
#'       
#' @param parameters list of intial parameter estimates
#' @param omega omega lower triangle (IIV)
#' @param err additive error
#' @param regimen `PKPDsim` regimen
#' @param model `PKPDsim` model
#' @param max_iter maximum number of iterations
#' @param data observed data
#' @param min_crit minimization criterion (default = 0.01)
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
  min_crit = 0.01,
  ...) {
  
  check_high_corr <- function(cov_mat, limit = 0.99) {
    tmp <- cov2cor(cov_mat)  
    if(any(tmp[lower.tri(tmp)] > limit)) {
      return(TRUE)
    } else {
      return(FALSE)
    }
  }

  ## initialize
  par_tmp   <- parameters
  omega_tmp <- omega
  err_tmp   <- err  # additive error
  cvs       <- diag(triangle_to_full(omega))
  ids       <- unique(data$id)
  message(paste0("Found ", length(ids), " subjects in dataset."))
  par_mn    <- unlist(parameters)
  all_pars  <- c(par_mn, diag(triangle_to_full(omega)))  
  
  ## loop iterations, only break when max_iter reached or stopping criterion reached
  for(j in 1:max_iter) {
    
    ## init and console output
    pars <- c()
    etas <- c()
    resid <- c()
    message(paste0("ITS iteration: ", j-1))
    message(" - Population means: ", paste(round(par_mn, 3), collapse = "  "))
    message(" - IIV %: ", paste(round(100 * sqrt(cvs), 1), collapse = "  "))
    message(paste0(" - Error: prop = ", err_tmp$prop, ", add = ", err_tmp$add))
    
    ## loop over patients
    pb <- txtProgressBar(min = 0, max = length(ids), initial = 0, style = 3)
    for(i in seq(ids)) {
      dv <- data[data$id == ids[i] && data$evid == 0,]
      dv <- data %>% filter(id == ids[i], evid == 0)
      tmp <- get_map_estimates(
        parameters = par_tmp,
        model = model,
        regimen = regimen,
        omega = omega_tmp,
        # weights = rep(1, length(dv[,1])),
        error = err_tmp,
        int_step_size = 1,
        data = dv,
        residuals = TRUE,
        include_error = TRUE,
        ...
      )
      setTxtProgressBar(pb, i)
      pars <- rbind(pars, cbind(t(as.numeric(unlist(tmp$parameters)))))
      etas <- rbind(etas, cbind(t(as.numeric(unlist(tmp$fit@coef)))))
      resid <- rbind(resid, cbind(t = tmp$t, pred = tmp$pred, ipred = tmp$ipred, ires = abs(tmp$ires), dv = tmp$dv))
    }
    
    ## update parameters to mean of estimates
    # par_mn <- apply(pars, 2, function(x) { exp(mean(log(x)) ) } )
    par_mn <- apply(pars, 2, mean)
    par_tmp <- parameters
    for(k in seq(names(parameters))) {
      par_tmp[[k]] <- par_mn[k]
    }
    
    ## update between-subject variability estimate
    omega_tmp1 <- cov(etas)
    if(!check_high_corr(omega_tmp1)) {
      omega_tmp <- omega_tmp1 # full_to_triangle(omega_tmp1)
    } else {
      message("High correlation in IIV matrix, resetting correlation to 0.")
      omega_tmp <- omega_tmp1 * diag(ncol(omega_tmp1))
    }
    cvs <- diag(omega_tmp1)

    ## update residual error estimate
    fit_comb <- glm(ires ~ ipred, data = data.frame(resid))
    intercept <- coef(fit_comb)[1]
    if(intercept < 0) {
      fit_prop <- glm(ires ~ ipred + 0 + offset(rep(intercept, length(resid[,1]))), data = data.frame(resid))
      err_tmp <- list(prop = coef(fit_prop)[1], add = 0)
    } else {
      err_tmp <- list(prop = coef(fit_comb)[2], add = coef(fit_comb)[1])
    }

    ## stopping criterion
    all_pars <- rbind(all_pars, c(par_mn, cvs))
    if(j > 1) {
      crit <- mean(abs(all_pars[j,] - all_pars[j-1,]) / all_pars[j-1,])
      message(paste0("\nStopping parameter: ", crit, " (min value ", min_crit, ")"))
      if(crit < min_crit) {
        message("Minimization criterion reached, stopping search.")
        return(list(parameters = all_pars, res = resid))
      }
    }
  }
  return(list(parameters = all_pars, res = resid))
} 
