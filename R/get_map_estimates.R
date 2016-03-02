#' Get MAP estimates
#'
#' @param model model
#' @param data data
#' @param parameters parameters
#' @param covariates covariates
#' @param fixed fix a specific parameters, supplied as vector of strings
#' @param weights vector of weights. Length of vector should be same as length of observation vector. If NULL, all weights are 1.
#' @param omega between subject variability, supplied as vector specifiying the lower triangle of the covariance matrix of random effects
#' @param error residual error, specified as list with arguments `add` and/or `prop` specifying the additive and proportional parts
#' @param regimen regimen
#' @param int_step_size integrator step size passed to PKPDsim
#' @param method optimization method, default L-BFGS-B
#' @param cols column names
#' @param residuals show residuals? This requires an additional simulation so will be slightly slower.
#' @param verbose show more output
#' @export
get_map_estimates <- function(
                      model = NULL,
                      data = NULL,
                      parameters = NULL,
                      covariates = NULL,
                      fixed = NULL,
                      weights = NULL,
                      omega = NULL,
                      error = list(prop = 0.1, add = 0.1, exp = 0),
                      regimen = NULL,
                      int_step_size = 0.1,
                      method = "L-BFGS-B",
                      type = "map",
                      cols = list(x="t", y="y"),
                      residuals = TRUE,
                      verbose = FALSE,
                      ...) {
  w_omega <- 1
  if(tolower(type) == "ls") {
    w_omega <- 0.001
  } # else either map or adaptive
  if(is.null(model) || is.null(data) || is.null(parameters) || is.null(omega) || is.null(regimen)) {
    stop("The 'model', 'data', 'omega', 'regimen', and 'parameters' arguments are required.")
  }
  if(!("function" %in% class(model))) {
    stop("The 'model' argument requires a function, e.g. a model defined using the new_ode_model() function from the PKPDsim package.")
  }
  if(!all(unlist(cols) %in% names(data))) {
    stop("Expected column names were not found in data. Please use 'cols' argument to specify column names for independent and dependent variable.")
  }
  if("PKPDsim" %in% class(data)) {
    if("comp" %in% names(data)) {
      data <- data[comp == "obs",]
      data <- data[!duplicated(data$t),]
      data$evid <- 0
    }
  }
  colnames(data) <- tolower(colnames(data))
  sig <- round(-log10(int_step_size))
  if(!("evid" %in% colnames(data))) {
    message("No 'evid' column in input data, assuming all rows are observations.")
  } else {
    data <- data[data$evid == 0,]
  }
  zero_offset <- NULL
  y_orig <- data$y
  if(any(data$t < min(regimen$dose_times))) { # protection against solving ODE from t < 0
    zero_offset <- rep(0, sum(data$t < min(regimen$dose_times)))
    filt <- data$t >= min(regimen$dose_times)
    if(!is.null(weights) && length(weights) == length(data$t) ) {
      weights <- weights[filt]
    }
    data <- data[filt,]
  }
  t_obs <- data$t
  if(!is.null(weights)) {
    if(length(weights) != length(t_obs)) {
      stop("Vector of weights of different size than observation vector!")
    }
  } else {
    weights <- 1
  }
  if(sum(unlist(error)) == 0) {
    stop("No residual error model specified, or residual error is 0.")
  }
  if(!is.null(attr(model, "cpp")) && attr(model, "cpp")) {
    ll_func <- function(
      data,
      eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9, eta10, eta11, eta12, # unfortunately seems no other way to do this...
      parameters,
      covariates,
      regimen = regimen,
      omega_full = omega_full,
      error = error,
      model,
      t_obs,
      sig,
      w_omega,
      covs) {
        par <- parameters
        p <- as.list(match.call())
        for(i in seq(names(par))) {
          par[[i]] <- par[[i]] * exp(p[[(paste0("eta", i))]])
        }
        suppressMessages({
          sim <- sim_ode(ode = model,
                         parameters = par,
                         covariates = covariates,
                         n_ind = 1,
                         int_step_size = int_step_size,
                         regimen = regimen,
                         t_obs = t_obs,
                         only_obs = TRUE,
                         ...)
        })
        ipred <- sim[!duplicated(sim$t),]$y
        y <- data$y
        res_sd <- sqrt(error$prop^2*ipred^2 + error$add^2)
        et <- mget(objects()[grep("eta", objects())])
        et <- as.numeric(as.character(et[et != ""]))
        omega_full <- omega_full[1:length(et), 1:length(et)]
        ofv <-   c(mvtnorm::dmvnorm(et, mean=rep(0, length(et)),
                                    sigma=omega_full[1:length(et), 1:length(et)],
                                    log=TRUE) * w_omega,
                   dnorm(y - ipred, mean = 0, sd = res_sd, log=TRUE) * weights)
        if(verbose) { print(ofv) }
        return(-sum(ofv))
      }
  } else {
    ll_func <- function(
      data,
      eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9, eta10, eta11, eta12, # unfortunately seems no other way to do this...
      parameters,
      covariates,
      regimen = regimen,
      omega_full = omega_full,
      error = error,
      model,
      t_obs,
      w_omega,
      sig,
      covs) {
        par <- parameters
        p <- as.list(match.call())
        for(i in seq(names(par))) {
          par[[i]] <- par[[i]] * exp(p[[(paste0("eta", i))]])
        }
        ipred <- model(t = data[data$evid == 0,]$t,
                       parameters = par)
        y <- data$y
        res_sd <- sqrt(error$prop^2*ipred^2 + error$add^2)
        et <- mget(objects()[grep("eta", objects())])
        et <- as.numeric(as.character(et[et != ""]))
        omega_full <- omega_full[1:length(et), 1:length(et)]
        ofv <-   c(mvtnorm::dmvnorm(et, mean=rep(0, length(et)),
                                    sigma=omega_full[1:length(et), 1:length(et)],
                                    log=TRUE) * w_omega,
                   dnorm(y - ipred, mean = 0, sd = res_sd, log=TRUE))
        if(verbose) { print(ofv) }
        return(-sum(ofv))
      }
  }
  eta <- list()
  for(i in seq(parameters)) {
    eta[[paste0("eta", i)]] <- 0
  }
  ## check if fixed parameter actually in parameter list
  fixed <- names(parameters)[names(parameters) %in% fixed]
  if(length(fixed) == 0) {
    fixed <- NULL
  }
  ## fix etas
  if (!is.null(fixed)) {
    id_fix <- match(fixed, names(parameters))
    fix <- list()
    for(i in 1:length(id_fix)) {
      id <- names(eta)[id_fix[i]]
      fix[[id]] <- 0
    }
  } else {
    fix <- NULL
  }
  omega_full <- diag(length(names(parameters))) # dummy om matrix
  om_nonfixed <- triangle_to_full(omega)
  if(nrow(om_nonfixed) < (length(parameters) - length(fix))) {
    stop("Provided omega matrix is smaller than expected based on the number of model parameters. Either fix some parameters or increase the size of the omega matrix.")
  }
  omega_full[1:nrow(om_nonfixed), 1:ncol(om_nonfixed)] <- om_nonfixed
  fit <- bbmle::mle2(ll_func,
              start = eta,
              method = method,
              data = list(data = data,
                          parameters = parameters,
                          covariates = covariates,
                          regimen = regimen,
                          model = model,
                          omega_full = omega_full,
                          error = error,
                          t_obs = t_obs,
                          w_omega = w_omega,
                          sig = sig,
                          covs = NULL),
              fixed = fix)
  cf <- bbmle::coef(fit)
  if(tolower(type) == "adaptive") { # not fully implemented yet.
    fit_ls <- bbmle::mle2(ll_func,
                          start = eta,
                          method = method,
                          data = list(data = data,
                                      parameters = parameters,
                                      covariates = covariates,
                                      regimen = regimen,
                                      model = model,
                                      omega_full = omega_full,
                                      error = error,
                                      t_obs = t_obs,
                                      w_omega = 0.001,
                                      sig = sig,
                                      covs = NULL),
                          fixed = fix)
    cf_ls <- bbmle::coef(fit_ls)
    if((cf_ls[1] / cf[1]) > 1.25) {
      cf <- cf_ls
    }
  }
  par <- parameters
  for(i in seq(names(par))) {
    par[[i]] <- as.numeric(as.numeric(par[[i]]) * exp(as.numeric(cf[i])))
  }
  if(residuals) {
    suppressMessages({
      sim_ipred <- sim_ode(ode = model,
                     parameters = par,
                     covariates = covariates,
                     n_ind = 1,
                     int_step_size = int_step_size,
                     regimen = regimen,
                     t_obs = t_obs,
                     only_obs = TRUE,
                     ...)
    })
    suppressMessages({
      sim_pred <- sim_ode(ode = model,
                     parameters = parameters,
                     covariates = covariates,
                     n_ind = 1,
                     int_step_size = int_step_size,
                     regimen = regimen,
                     t_obs = t_obs,
                     only_obs = TRUE,
                     ...)
    })
    ipred <- sim_ipred[!duplicated(sim_ipred$t),]$y
    pred <- sim_pred[!duplicated(sim_pred$t),]$y
    w_ipred <- sqrt(error$prop^2 * ipred^2 + error$add^2)
    w_pred <- sqrt(error$prop^2 * pred^2 + error$add^2)
    y <- data$y
    res <- (y - pred)
    wres <- res / w_pred # not same as wres in NONMEM!
    ires <- (y - ipred)
    iwres <- ires / w_ipred
  }
  obj <- list(fit = fit, parameters = par)
  if(residuals) {
    obj$res <- c(zero_offset, res)
    obj$wres <- c(zero_offset, wres)
    obj$ires <- c(zero_offset, ires)
    obj$iwres <- c(zero_offset, iwres)
    obj$ipred <- c(zero_offset, ipred)
    obj$pred <- c(zero_offset, pred)
    obj$dv <- y_orig
  }
  class(obj) <- c(class(obj), "map_estimates")
  return(obj)
}
