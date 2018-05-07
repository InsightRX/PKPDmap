#' Get MAP estimates
#'
#' @param model model, created using `PKPDsim::new_ode_model()`
#' @param data data data.frame with columns `t` and `y` (and possibly evid)
#' @param parameters list of parameters
#' @param covariates list of covariates, each one created using `PKPDsim::new_coviarate()`
#' @param fixed fix a specific parameters, supplied as vector of strings
#' @param as_eta vector of parameters that are estimates as eta (e.g. IOV)
#' @param weights vector of weights for error. Length of vector should be same as length of observation vector. If NULL (default), all weights are equal. Used in both MAP and NP methods. Note that `weights` argument will also affect residuals (residuals will be scaled too).
#' @param omega between subject variability, supplied as vector specifiying the lower triangle of the covariance matrix of random effects
#' @param weight_prior weighting of priors in relationship to observed data, default = 1
#' @param iov_bins bins for inter-occasion variability. Passed unchanged to PKPDsim.
#' @param error residual error, specified as list with arguments `add` and/or `prop` specifying the additive and proportional parts
#' @param ltbs log-transform both sides? (`NULL` by default, meaning that it will be picked up from the PKPDsim model. Can be overridden with `TRUE`). Note: `error` should commonly only have additive part.
#' @param censoring label for column specifying censoring. If value in dataset in this column is < 0 then censoring is assumed <LLOQ. If > 0 then  >ULOQ.
#' @param mixture specify mixture model. Currently for single parameter only. Overwrites regular parameter, if specified. Specify e.g. as: `mixture = list("CL" = c(5, 9))`
#' @param include_omega TRUE
#' @param include_error TRUE
#' @param regimen regimen
#' @param A_init initial state vector
#' @param int_step_size integrator step size passed to PKPDsim
#' @param optimizer optimization library to use, default is `optim`
#' @param method optimization method, default `BFGS`
#' @param control list of options passed to `optim()` function
#' @param type estimation type, options are `map`, `ls`, and `np_hybrid`
#' @param np_settings list with settings for non-parametric estimation (if selected), containing any of the following: `error`, `grid_span`, grid_size`, `grid_exponential`
#' @param cols column names
#' @param residuals show residuals? This requires an additional simulation so will be slightly slower.
#' @param output_include passed to PKPDsim::sim_ode(), returns covariates and parmeter values over time in return object. Only invoked if `residuals` option is `TRUE`.
#' @param skip_hessian skip calculation of Hessian (in `bbmle::mle2()`)
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
                      iov_bins = NULL,
                      error = NULL,
                      ltbs = NULL,
                      censoring = NULL,
                      mixture = NULL,
                      weights = NULL,
                      include_omega = TRUE,
                      include_error = TRUE,
                      regimen = NULL,
                      t_init = 0,
                      int_step_size = 0.1,
                      optimizer = "optim",
                      method = "BFGS",
                      control = list(reltol = 1e-5),
                      type = "map",
                      np_settings = list(),
                      cols = list(x = "t", y = "y"),
                      residuals = TRUE,
                      verbose = FALSE,
                      A_init = NULL,
                      skip_hessian = FALSE,
                      output_include = list(covariates = FALSE, parameters = FALSE),
                      ...) {

  if(optimizer == "optimx") require("optimx")

  ## Handle weighting of priors, allow for some presets but can
  ## also be set manually using `weight_prior`
  if(is.null(weight_prior) || is.na(weight_prior)) {
    weight_prior <- 1
  }
  weight_prior <- weight_prior^2 # 3x larger IIV on SD scale
  calc_ofv <- calc_ofv_map
  if(weight_prior == 0) {
    calc_ofv <- calc_ofv_ls
  }
  if(tolower(type) %in% c("map", "pls")) {
    if(is.null(model) || is.null(data) || is.null(parameters) || is.null(omega) || is.null(regimen)) {
      stop("The 'model', 'data', 'omega', 'regimen', and 'parameters' arguments are required.")
    }
  }
  if(tolower(type) == "pls") {
    weight_prior <- 0.001
    ## RK: Empirically determined to be a good penalty.
    ##     In principle, PLS is just MAP with very flat priors
  }
  if(tolower(type) %in% c("ls")) {
    if(is.null(model) || is.null(data) || is.null(regimen)) {
      stop("The 'model', 'data', and 'parameters' arguments are required.")
    }
    calc_ofv <- calc_ofv_ls
    error <- list(prop = 0, add = 1)
  }
  if(!is.null(error)) { ## safety checks
    if(is.null(error$prop)) error$prop <- 0
    if(is.null(error$add)) error$add <- 0
    if(is.null(error$exp)) error$exp <- 0
  } else {
    warning("No residual error specified, using: 10% proportional + 0.1 additive error.")
    list(prop = 0.1, add = 0.1, exp = 0)
  }
  if(!is.null(censoring) && class(censoring) != "character") {
    stop("Censoring argument requires label specifying column in dataset with censoring info.")
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
      # data <- data[!duplicated(data$t),]
      data$evid <- 0
    }
  }
  colnames(data) <- tolower(colnames(data))
  sig <- round(-log10(int_step_size))
  if("evid" %in% colnames(data)) {
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
  if(any(duplicated(t_obs))) message("Duplicate times were detected in data. Estimation will proceed but please check that data is correct. For putting more weight on certain measurements, please use the `weights` argument.")
  if(!is.null(weights)) {
    if(length(weights) != length(t_obs)) {
      stop("Vector of weights of different size than observation vector!")
    }
  } else {
    weights <- rep(1, length(data$y))
  }
  if(sum(unlist(error)) == 0) {
    stop("No residual error model specified, or residual error is 0.")
  }

  ## Log-transform both sides?
  if(is.null(ltbs)) {
    ltbs <- FALSE
    if(!is.null(attr(model, "ltbs")) && attr(model, "ltbs")) ltbs <- TRUE
  }
  if(ltbs) {
    transf <- function(x) log(x)
  } else {
    transf <- function(x) x
  }

  #################################################
  ## Likelihood function using PKPDsim
  #################################################
  ll_func_PKPDsim <- function(
    data,
    sim_object,
    parameters,
    t_obs,
    # unfortunately seems no other way to do this using the optim package...
    eta01, eta02, eta03, eta04, eta05, eta06, eta07, eta08, eta09, eta10,
    eta11, eta12, eta13, eta14, eta15, eta16, eta17, eta18, eta19, eta20,
    eta21, eta22, eta23, eta24,
    error = error,
    model,
    sig,
    weight_prior,
    as_eta,
    censoring_idx,
    censoring_label,
    iov_bins,
    ...) {
    par <- parameters
    p <- as.list(match.call())
    for(i in seq(nonfixed)) {
      key <- nonfixed[i]
      if(key %in% p$as_eta) {
        par[[key]] <- p[[(paste0("eta", sprintf("%02d", i)))]]
      } else {
        par[[key]] <- par[[key]] * exp(p[[(paste0("eta", sprintf("%02d", i)))]])
      }
    }
    sim_object$p <- par
    ipred <- transf(PKPDsim::sim_core(
      sim_object,
      ode = model,
      duplicate_t_obs = TRUE,
      t_init = t_init)$y)
    dv <- transf(data$y)
    ofv_cens <- NULL
    if(!is.null(censoring_idx)) {
      cens <- data[[censoring_label]][censoring_idx] # if <0 then LLOQ, if >0 ULOQ
      cens <- ifelse(cens > 0, -1, 1)
      ipred_cens <- ipred[censoring_idx]
      ipred <- ipred[!censoring_idx]
      dv_cens <- dv[censoring_idx]
      dv <- dv[!censoring_idx]
      weights_cens <- weights[censoring_idx]
      weights <- weights[!censoring_idx]
      res_sd_cens <- sqrt(error$prop^2*ipred_cens^2 + error$add^2)
      ofv_cens <- stats::pnorm((dv_cens - ipred_cens) * cens, 0, res_sd_cens, log=TRUE) * weights_cens
    }
    res_sd <- sqrt(error$prop^2*ipred^2 + error$add^2)
    et <- mget(objects()[grep("^eta", objects())])
    et <- as.numeric(as.character(et[et != ""]))
    omega_full <- as.matrix(omega_full)[1:length(et), 1:length(et)]
    ofv <- calc_ofv(
      eta = et,
      omega = omega_full,
      dv = dv,
      ipred = ipred,
      res_sd = res_sd,
      weights = weights,
      weight_prior = weight_prior,
      include_omega = include_omega,
      include_error = include_error)
    ofv <- c(ofv, ofv_cens)
    if(verbose) {
      cat("-------------------------------------------------------------\n")
      cat(paste0("Eta\t: [", paste(signif(et,5), collapse=", "),"]\n"))
      cat(paste0("y_hat\t: [", paste(ipred, collapse=", "),"]\n"))
      cat(paste0("P(y)\t: [", paste(signif(exp(ofv[-1]),5), collapse=", "),"]\n"))
      cat(paste0("OFV\t: [", paste(signif(-2*sum(ofv),5), collapse=", "), "]\n"))
    }
    return(-2 * sum(ofv))
  }

  ll_func <- ll_func_PKPDsim
  if(is.null(attr(model, "cpp")) || !attr(model, "cpp")) {
    ll_func <- ll_func_generic
  }

  ## check if fixed parameter actually in parameter list
  if(length(intersect(fixed, names(parameters))) != length(fixed)) {
    warning("Warning: not all fixed parameters were found in parameter set!\n")
  }
  fixed <- names(parameters)[names(parameters) %in% fixed]
  if(length(fixed) == 0) {
    fixed <- NULL
  }

  ## fix etas
  nonfixed <- names(parameters)[is.na(match(names(parameters), fixed))]
  n_nonfix <- length(nonfixed)
  eta <- list()
  for(i in seq(nonfixed)) {
    eta[[paste0("eta", sprintf("%02d", i))]] <- 0
  }
  fix <- NULL

#  omega_full <- diag(length(names(parameters))) # dummy om matrix
  if(class(omega) == "matrix") {
    omega_full <- omega # dummy om matrix
  } else {
    omega_full <- triangle_to_full(omega) # dummy om matrix
  }
  om_nonfixed <- triangle_to_full(omega)
  if(nrow(om_nonfixed) < (length(parameters) - length(fixed))) {
    msg <- "Provided omega matrix is smaller than expected based on the number of model parameters. Either fix some parameters or increase the size of the omega matrix.\n"
    msg <- c(msg,
      paste0("Non-fixed omegas: ", paste(om_nonfixed, collapse=", "), "\n"),
      paste0("Parameters: ", paste(parameters, collapse=", "), "\n"),
      paste0("Fixed: ", paste(fix, collapse=","), "\n"))
    stop(msg)
  }

  ## check if censoring code needs to be used
  censoring_idx <- NULL
  if(!is.null(censoring)) {
    if(any(data[[tolower(censoring)]] != 0)) {
      censoring_idx <- data[[tolower(censoring)]] != 0
      if(verbose) message("One or more values in data are censored, including censoring in likelihood.")
    } else {
      if(verbose) message("Warning: censoring specified, but no censored values in data.")
    }
  }

  #################################################
  ## create simulation design up-front:
  #################################################
  suppressMessages({
    sim_object <- PKPDsim::sim(ode = model,
                               parameters = parameters,
                               covariates = covariates,
                               n_ind = 1,
                               int_step_size = int_step_size,
                               regimen = regimen,
                               t_obs = t_obs,
                               checks = FALSE,
                               only_obs = TRUE,
                               A_init = A_init,
                               fixed = fixed,
                               t_max = tail(t_obs, 1) + t_init + 1,
                               iov_bins = iov_bins,
                               return_design = TRUE,
                               t_init = t_init,
                               ...)
  })

  mixture_obj <- NULL
  if(!is.null(mixture)) {
    if((class(mixture) != "list") || length(names(mixture)) != 1) {
      stop("`mixture` argument needs to be a list containing elements `parameters` and `values`.")
    }
    mix_par <- names(mixture)[1]
    mix_par_values <- mixture[[mix_par]]
    fits <- list()
    par_mix <- parameters
    ofvs <- c()
    for(i in seq(mix_par_values)) {
      par_mix[[mix_par]] <- mix_par_values[i]
      fits[[i]] <- bbmle::mle2(ll_func,
                         start = eta,
                         method = method,
                         optimizer = optimizer,
                         control = control,
                         skip.hessian = skip_hessian,
                         data = list(data = data,
                                     sim_object = sim_object,
                                     parameters = par_mix,
                                     t_obs = t_obs,
                                     model = model,
                                     error = error,
                                     weight_prior = weight_prior,
                                     sig = sig,
                                     as_eta = as_eta,
                                     censoring_idx = censoring_idx,
                                     censoring_label = censoring,
                                     iov_bins = iov_bins),
                         fixed = fix)
      ofvs <- c(ofvs, bbmle::logLik(fits[[i]]))
    }
    like <- exp(ofvs)
    prob <- like / sum(like)
    fit <- fits[[match(max(prob), prob)]]
    mixture_obj <- list(
      parameter = mix_par,
      values = mix_par_values,
      selected = mix_par_values[match(max(prob), prob)],
      probabilities = prob
    )
  } else {
    output <- tryCatch({
      fit <- bbmle::mle2(ll_func,
                       start = eta,
                       method = method,
                       optimizer = optimizer,
                       control = control,
                       skip.hessian = skip_hessian,
                       data = list(data = data,
                                   sim_object = sim_object,
                                   parameters = parameters,
                                   t_obs = t_obs,
                                   model = model,
                                   error = error,
                                   weight_prior = weight_prior,
                                   sig = sig,
                                   as_eta = as_eta,
                                   censoring_idx = censoring_idx,
                                   censoring_label = censoring,
                                   iov_bins = iov_bins),
                       fixed = fix)
    }, error = function(e) {
       return(e)
    })
    if("error" %in% class(output)) return(output)
  }
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
  obj <- list(fit = fit, mixture = mixture_obj)

  #################################################
  ## Non-parametric estimation
  #################################################
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

  #################################################
  ## Add g.o.f. info
  #################################################
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
                           output_include = output_include,
                           t_init = t_init,
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
                          t_init = t_init,
                          ...)
    })
    ipred <- sim_ipred$y
    pred <- sim_pred$y
    w_ipred <- sqrt(error$prop^2 * transf(ipred)^2 + error$add^2)
    w_pred <- sqrt(error$prop^2 * transf(pred)^2 + error$add^2)
    y <- data$y
    prob <- list(par = c(mvtnorm::pmvnorm(cf, mean=rep(0, length(cf)),
                         sigma = omega_full[1:length(cf), 1:length(cf)])),
                 data = stats::pnorm(transf(y) - transf(ipred), mean = 0, sd = w_ipred))
    res <- (transf(y) - transf(pred))
    wres <- (res / w_pred) * weights
    cwres <- res / sqrt(abs(cov(transf(pred), transf(y)))) * weights
    # Note: in NONMEM CWRES is on the population level, so can't really compare. NONMEM calls this CIWRES, it seems.
    ires <- (transf(y) - transf(ipred))
    iwres <- (ires / w_ipred)
    iwres_weighted <- iwres * weights
    obj$prob <- prob
    if(length(w_ipred) > 1) {
      obj$mahalanobis <- stats::mahalanobis(transf(y), transf(ipred), cov = diag(w_ipred^2))
    } else {
      obj$mahalanobis <- stats::mahalanobis(transf(y), transf(ipred), cov = w_ipred^2)
    }
    obj$res <- c(zero_offset, res)
    obj$wres <- c(zero_offset, cwres)
    obj$cwres <- c(zero_offset, cwres)
    obj$ires <- c(zero_offset, ires)
    obj$iwres <- c(zero_offset, iwres)
    obj$iwres_weighted <- c(zero_offset, iwres_weighted)
    obj$ipred <- c(zero_offset, ipred)
    obj$pred <- c(zero_offset, pred)
    obj$weights <- c(zero_offset, weights)
    obj$dv <- y_orig
    if(output_include$covariates && !is.null(covariates)) {
      obj$covariates_time <- sim_ipred[!duplicated(sim_ipred$t), names(covariates)]
    }
    if(output_include$parameters) {
      obj$parameters_time <- sim_ipred[!duplicated(sim_ipred$t), names(parameters)]
    }
  }
  obj$vcov_full <- fit@vcov
  obj$vcov <- fit@vcov[t(!upper.tri(fit@vcov))]
  class(obj) <- c(class(obj), "map_estimates")
  return(obj)
}
