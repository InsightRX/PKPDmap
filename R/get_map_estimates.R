#' Get MAP estimates
#'
#' @param model model, created using `PKPDsim::new_ode_model()`
#' @param data data data.frame with columns `t` and `y` (and possibly evid)
#' @param parameters list of parameters
#' @param covariates list of covariates, each one created using `PKPDsim::new_coviarate()`
#' @param fixed fix a specific parameters, supplied as vector of strings
#' @param as_eta vector of parameters that are estimates as eta (e.g. IOV)
#' @param weights vector of weights for error. Length of vector should be same as length of observation vector. If NULL (default), all weights are equal. Used in both MAP and NP methods.
#' @param omega between subject variability, supplied as vector specifiying the lower triangle of the covariance matrix of random effects
#' @param weight_prior weighting of priors in relationship to observed data, default = 1
#' @param error residual error, specified as list with arguments `add` and/or `prop` specifying the additive and proportional parts
#' @param include_omega TRUE
#' @param include_error TRUE
#' @param regimen regimen
#' @param int_step_size integrator step size passed to PKPDsim
#' @param method optimization method, default BFGS
#' @param type estimation type, options are `map`, `ls`, and `np_hybrid`
#' @param np_settings list with settings for non-parametric estimation (if selected), containing any of the following: `error`, `grid_span`, grid_size`, `grid_exponential`
#' @param cols column names
#' @param residuals show residuals? This requires an additional simulation so will be slightly slower.
#' @param verbose show more output
#' @param ... parameters passed on to `sim_ode()` function
#' @export
get_map_estimates <- function(
                      model = NULL,
                      data = NULL,
                      parameters = NULL,
                      covariates = NULL,
                      fixed = c(),
                      as_eta = c(),
                      omega = NULL,
                      weight_prior = 1,
                      error = list(prop = 0.1, add = 0.1, exp = 0),
                      weights = NULL,
                      include_omega = TRUE,
                      include_error = TRUE,
                      regimen = NULL,
                      int_step_size = 0.1,
                      method = "BFGS",
                      type = "map",
                      np_settings = list(),
                      cols = list(x = "t", y = "y"),
                      residuals = TRUE,
                      verbose = FALSE,
                      A_init = NULL,
                      ...) {

  ## Handle weighting of priors, allow for some presets but can
  ## also be set manually using `weight_prior`
  if(is.null(weight_prior) || is.na(weight_prior)) {
    weight_prior <- 1
  }
  weight_prior <- weight_prior^2 # 3x larger IIV on SD scale
  if(tolower(type) == "ls") {
    weight_prior <- 0.001
  }
  if(!is.null(error)) { ## safety checks
    if(is.null(error$prop)) error$prop <- 0
    if(is.null(error$add)) error$add <- 0
    if(is.null(error$exp)) error$exp <- 0
  }

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
      data <- data[data$comp == "obs",]
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

  ## Likelihood function for PKPDsim
  ll_func_PKPDsim <- function(
    data,
    # unfortunately seems no other way to do this...
    eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9, eta10,
    eta11, eta12, eta13, eta14, eta15, eta16, eta17, eta18, eta19, eta20,
    eta21, eta22, eta23, eta24,
    parameters,
    covariates = NULL,
    covariate_names = NULL,
    regimen = regimen,
    omega_full = omega_full,
    error = error,
    model,
    t_obs,
    sig,
    int_step_size = 0.1,
    weight_prior,
    as_eta,
    ...) {
    par <- parameters
    if(!is.null(covariates)) { # not properly passed through bbmle it seems
      for (i in seq(covariate_names)) {
        names(covariates)[i] <- covariate_names[i]
      }
    }
    p <- as.list(match.call())
    for(i in seq(names(par))) {
      key <- names(par)[i]
      if(key %in% as_eta) {
        par[[key]] <- par[[key]] + p[[(paste0("eta", i))]]
      } else {
        par[[key]] <- par[[key]] * exp(p[[(paste0("eta", i))]])
      }
    }
    suppressMessages({
      sim <- PKPDsim::sim_ode(ode = model,
                              parameters = par,
                              covariates = covariates,
                              n_ind = 1,
                              int_step_size = int_step_size,
                              regimen = regimen,
                              t_obs = t_obs,
                              checks = FALSE,
                              only_obs = TRUE,
                              A_init = A_init,
                              ...)
    })
    ipred <- sim[!duplicated(sim$t),]$y
    res_sd <- sqrt(error$prop^2*ipred^2 + error$add^2)
    et <- mget(objects()[grep("^eta", objects())])
    et <- as.numeric(as.character(et[et != ""]))
    omega_full <- omega_full[1:length(et), 1:length(et)]
    ofv <-   c(mvtnorm::dmvnorm(et, mean=rep(0, length(et)),
                                sigma = omega_full[1:length(et),
                                                   1:length(et)] * 1/weight_prior,
                                log=TRUE) * include_omega,
               stats::dnorm((data$y - ipred) * include_error, mean = 0, sd = res_sd, log=TRUE) * weights)
    if(verbose) {
      print(ofv)
    }
    return(-sum(ofv))
  }

  ## generic likelihood function
  ll_func_generic <- function(
    data,
    # unfortunately seems no other way to do this...
    eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9, eta10,
    eta11, eta12, eta13, eta14, eta15, eta16, eta17, eta18, eta19, eta20,
    eta21, eta22, eta23, eta24,
    parameters,
    covariates = NULL,
    covariate_names = NULL,
    regimen = regimen,
    omega_full = omega_full,
    error = error,
    model,
    t_obs,
    sig,
    weight_prior,
    ...) {
    par <- parameters
    p <- as.list(match.call())
    if(!is.null(covariates)) { # not properly passed through bbmle it seems
      for (i in seq(covariate_names)) {
        names(covariates)[i] <- covariate_names[i]
      }
    }
    for(i in seq(names(par))) {
      par[[i]] <- par[[i]] * exp(p[[(paste0("eta", i))]])
    }
    ipred <- model(t = data[data$evid == 0,]$t,
                   parameters = par)
    y <- data$y
    res_sd <- sqrt(error$prop^2*ipred^2 + error$add^2)
    et <- mget(objects()[grep("^eta", objects())])
    et <- as.numeric(as.character(et[et != ""]))
    omega_full <- omega_full[1:length(et), 1:length(et)]
    ofv <-   c(mvtnorm::dmvnorm(et, mean=rep(0, length(et)),
                                sigma = omega_full[1:length(et), 1:length(et)],
                                log=TRUE) * weight_prior,
               stats::dnorm(y - ipred, mean = 0, sd = res_sd, log=TRUE) * weights)
    if(verbose) {
      print(ofv)
    }
    return(-sum(ofv))
  }

  ll_func <- ll_func_PKPDsim
  if(is.null(attr(model, "cpp")) || !attr(model, "cpp")) {
    ll_func <- ll_func_generic
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
  nonfixed <- names(parameters)[is.na(match(names(parameters), fixed))]
  n_nonfix <- length(nonfixed)
  if (!is.null(fixed)) {
    if(n_nonfix == 0) {
      stop("Nothing to estimate!")
    }
    id_fix <- match(fixed, names(parameters))
    fix <- list()
    for(i in (n_nonfix+1):length(parameters)) {
      id <- names(eta)[i]
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
  omega_full[1:n_nonfix, 1:n_nonfix] <- om_nonfixed[1:n_nonfix, 1:n_nonfix]
  fit <- bbmle::mle2(ll_func,
              start = eta,
              method = method,
              data = list(data = data,
                          parameters = parameters,
                          covariates = covariates,
                          covariate_names = names(covariates),
                          A_init = A_init,
                          regimen = regimen,
                          model = model,
                          omega_full = omega_full,
                          error = error,
                          t_obs = t_obs,
                          weight_prior = weight_prior,
                          sig = sig,
                          as_eta = as_eta,
                          int_step_size = int_step_size),
              fixed = fix)
  cf <- bbmle::coef(fit)
  par <- parameters
  for(i in seq(nonfixed)) {
    key <- nonfixed[i]
    if(key %in% as_eta) {
      par[[key]] <- as.numeric(cf[i])
    } else {
      par[[key]] <- as.numeric(as.numeric(par[[key]]) * exp(as.numeric(cf[i])))
    }
  }
  obj <- list(fit = fit)
  if(type == "np_hybrid") {
    obj$parameters_map <- par ## keep MAP estimates
    np_settings <- replace_list_elements(np_settings_default, np_settings)
    if(is.null(np_settings$error)) { # if no specific error magnitude is specified for NP, just use same as used for MAP
      np_settings$error <- error
    }
    if(np_settings$grid_adaptive) { # do a first pass with a broad grid
      pars_grid <- create_grid_around_parameters(
        parameters = par,
        span = np_settings$grid_span_adaptive,
        exponential = np_settings$grid_exponential_adaptive,
        grid_size = np_settings$grid_size_adaptive)
      np <- get_np_estimates(parameter_grid = pars_grid,
                             error = np_settings$error,
                             model = model,
                             regimen = regimen,
                             data = data$y,
                             t_obs = data$t,
                             covariates = covariates,
                             weights = weights)
      # take the estimates with highest probability as starting point for next grid
      tmp <- np$prob[order(-np$prob$like),][1,]
      par <- as.list(tmp[1:length(np$parameters)])
    }
    pars_grid <- create_grid_around_parameters(
      parameters = par,
      span = np_settings$grid_span,
      exponential = np_settings$grid_exponential,
      grid_size = np_settings$grid_size)
    np <- get_np_estimates(parameter_grid = pars_grid,
                           error = np_settings$error,
                           model = model,
                           regimen = regimen,
                           data = data$y,
                           t_obs = data$t,
                           covariates = covariates,
                           weights = weights)
    for(i in 1:length(par)) {
      par[[i]] <- np$parameters[[i]]
    }
    obj$np <- list(prob = np$prob)
  }
  obj$parameters <- par
  if(residuals) {
    suppressMessages({
      sim_ipred <- PKPDsim::sim_ode(ode = model,
                           parameters = par,
                           covariates = covariates,
                           n_ind = 1,
                           int_step_size = int_step_size,
                           regimen = regimen,
                           t_obs = t_obs,
                           only_obs = TRUE,
                           checks = FALSE,
                           A_init = A_init,
                           ...)
    })
    suppressMessages({
      sim_pred <- PKPDsim::sim_ode(ode = model,
                          parameters = parameters,
                          covariates = covariates,
                          n_ind = 1,
                          int_step_size = int_step_size,
                          regimen = regimen,
                          t_obs = t_obs,
                          only_obs = TRUE,
                          checks = FALSE,
                          A_init = A_init,
                          ...)
    })
    ipred <- sim_ipred[!duplicated(sim_ipred$t),]$y
    pred <- sim_pred[!duplicated(sim_pred$t),]$y
    w_ipred <- sqrt(error$prop^2 * ipred^2 + error$add^2)
    w_pred <- sqrt(error$prop^2 * pred^2 + error$add^2)
    y <- data$y
    prob <- list(par = c(mvtnorm::pmvnorm(cf, mean=rep(0, length(cf)),
                         sigma = omega_full[1:length(cf), 1:length(cf)])),
                 data = stats::pnorm(y - ipred, mean = 0, sd = w_ipred))
    res <- (y - pred)
    wres <- res / w_pred
    cwres <- res / sqrt(cov(pred, y)) # Note: in NONMEM this is on the population level, so can't really compare
    ires <- (y - ipred)
    iwres <- ires / w_ipred
    obj$prob <- prob
    if(length(w_ipred) > 1) {
      obj$mahalanobis <- stats::mahalanobis(y, ipred, cov = diag(w_ipred^2))
    } else {
      obj$mahalanobis <- stats::mahalanobis(y, ipred, cov = w_ipred^2)
    }
    obj$res <- c(zero_offset, res)
    obj$wres <- c(zero_offset, wres)
    obj$wres <- c(zero_offset, cwres)
    obj$ires <- c(zero_offset, ires)
    obj$iwres <- c(zero_offset, iwres)
    obj$ipred <- c(zero_offset, ipred)
    obj$pred <- c(zero_offset, pred)
    obj$dv <- y_orig
  }
  obj$vcov_full <- fit@vcov
  obj$vcov <- fit@vcov[t(!upper.tri(fit@vcov))]
  class(obj) <- c(class(obj), "map_estimates")
  return(obj)
}
