#' Get MAP estimates
#'
#' @param model model, created using `PKPDsim::new_ode_model()`
#' @param data data data.frame with columns `t` and `y` (and possibly evid)
#' @param parameters list of parameters
#' @param covariates list of covariates, each one created using `PKPDsim::new_coviarate()`
#' @param fixed fix a specific parameters, supplied as vector of strings
#' @param as_eta vector of parameters that are estimates as eta (e.g. IOV)
#' @param weights vector of weights for error. Length of vector should be same 
#' as length of observation vector. If NULL (default), all weights are equal. 
#' Used in both MAP and NP methods. Note that `weights` argument will also 
#' affect residuals (residuals will be scaled too).
#' @param omega between subject variability, supplied as vector specifiying the
#'  lower triangle of the covariance matrix of random effects
#' @param weight_prior weighting of priors in relationship to observed data, 
#' default = 1
#' @param iov_bins bins for inter-occasion variability. Passed unchanged to 
#' PKPDsim.
#' @param error residual error, specified as list with arguments `add` and/or 
#' `prop` specifying the additive and proportional parts
#' @param ltbs log-transform both sides? (`NULL` by default, meaning that it 
#' will be picked up from the PKPDsim model. Can be overridden with `TRUE`). 
#' Note: `error` should commonly only have additive part.
#' @param obs_type_label column name in `data` referring to observation type. 
#' Can be used for specification of different residual error models for 
#' differing observation types (e.g. venous and capillary or parent and 
#' metabolite), Residual error should then be specified as list of vectors, 
#' e.g. `list(prop = c(0.2, 0.1), add = c(1, 2))`.
#' @param censoring label for column specifying censoring. If value in dataset 
#' in this column is < 0 then censoring is assumed <LLOQ. If > 0 then  >ULOQ.
#' @param steady_state_analytic list object with settings for steady state MAP 
#' estimation.
#' @param include_omega TRUE
#' @param include_error TRUE
#' @param regimen regimen
#' @param t_init initialization time before first dose, default 0.
#' @param A_init initial state vector
#' @param int_step_size integrator step size passed to PKPDsim
#' @param ll_func likelihood function, default is `ll_func_PKPDsim` as included 
#' in this package.
#' @param optimizer optimization library to use, default is `optim`
#' @param method optimization method, default `BFGS`
#' @param control list of options passed to `optim()` function
#' @param allow_obs_before_dose allow observation before first dose?
#' @param type estimation type, options are `map`, `ls`, and `np_hybrid`
#' @param np_settings list with settings for non-parametric estimation (if 
#' selected), containing any of the following: `error`, `grid_span`, grid_size`,
#'  `grid_exponential`
#' @param cols column names
#' @param residuals show residuals? This requires an additional simulation so 
#' will be slightly slower.
#' @param output_include passed to PKPDsim::sim_ode(), returns covariates and
#'parameter values over time in return object. Only invoked if `residuals` 
#'option is `TRUE`.
#' @param skip_hessian skip calculation of Hessian
#' @param verbose show more output
#' @param ... parameters passed on to `sim_ode()` function
#' 
#' @examples
#' \dontrun{
#' ## define parameters
#' pk1 <- new_ode_model(
#'   code = "dAdt[1] = -(CL/V)*A[1]", 
#'   obs = list(scale = "V", cmt=1)
#' )
#' regimen <- new_regimen(
#'   amt = 100, interval = 12, times=c(0,24), 
#'   type="infusion"
#' )
#' parameters <- list("CL" = 3, "V" = 10) 
#' omega <- c(0.1, 
#'            0.02, 0.1)
#' 
#' ## simulate single individual
#' data <- sim_ode(
#'   ode = pk1, 
#'   parameters = list(CL = 6, V = 14), 
#'   regimen = regimen,
#'   omega = omega,
#'   t_obs = c(23, 25, 27, 32, 36)
#' )
#' 
#' ## Fit
#' get_map_estimates(
#'   model = pk1, 
#'   data = data, 
#'   parameters = parameters,
#'   omega = omega,
#'   regimen = regimen,
#'   error = list(prop = 0.1, add = 0.1)
#' )
#' }
#' 
#' @export
#' 
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
                      obs_type_label = NULL,
                      censoring = NULL,
                      weights = NULL,
                      steady_state_analytic = NULL,
                      include_omega = TRUE,
                      include_error = TRUE,
                      regimen = NULL,
                      t_init = 0,
                      int_step_size = 0.1,
                      ll_func = ll_func_PKPDsim,
                      optimizer = "optim",
                      method = "BFGS",
                      control = list(reltol = 1e-5),
                      allow_obs_before_dose = FALSE,
                      type = "map",
                      np_settings = list(),
                      cols = list(x = "t", y = "y"),
                      residuals = TRUE,
                      verbose = FALSE,
                      A_init = NULL,
                      skip_hessian = FALSE,
                      output_include = list(covariates = FALSE, parameters = FALSE),
                      ...) {

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
  if(is.null(obs_type_label)) {
    data$obs_type <- 1
  } else {
    data$obs_type <- data[[obs_type_label]]
  }
  if(!is.null(error)) { ## safety checks
    if(is.null(error$prop)) error$prop <- 0
    if(is.null(error$add)) error$add <- 0
    if(is.null(error$exp)) error$exp <- 0
  } else {
    warning("No residual error specified, using: 10% proportional + 0.1 additive error.")
    list(prop = 0.1, add = 0.1, exp = 0)
  }
  if(!is.null(censoring) && !inherits(censoring, "character")) {
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
  data_before_init <- data.frame()
  data <- data[order(data$t, data$obs_type),]
  y_orig <- data$y
  if(!allow_obs_before_dose) {
    filt_before_init <- data$t < min(regimen$dose_times)
    data_before_init <- data[filt_before_init,]
    if(any(filt_before_init)) { # protection against solving ODE from t < 0
      if(!is.null(weights) && length(weights) == length(data$t) ) {
        weights <- weights[!filt_before_init]
      }
      data <- data[!filt_before_init,]
    }
  }
  t_obs <- data$t
  if(any(duplicated(paste(t_obs, data$obs_type, sep = "_")))) {
    message("Duplicate times were detected in data. Estimation will proceed but please check that data is correct. For putting more weight on certain measurements, please use the `weights` argument.")
  }
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
  mixture_group <- 1
  suppressMessages({
    sim_object <- PKPDsim::sim(
      ode = model,
      parameters = parameters,
      covariates = covariates,
      n_ind = 1,
      int_step_size = int_step_size,
      regimen = regimen,
      t_obs = t_obs,
      obs_type = data$obs_type,
      checks = FALSE,
      only_obs = TRUE,
      A_init = A_init,
      fixed = fixed,
      mixture_group = mixture_group, # dummy value
      t_max = tail(t_obs, 1) + t_init + 1,
      iov_bins = iov_bins,
      return_design = TRUE,
      t_init = t_init,
      ...
    )
  })

  omega_full_est <- omega_full[1:n_nonfix, 1:n_nonfix]
  mixture_obj <- NULL
  if(!is.null(attr(model, "mixture"))) {
    mixture <- attr(model, "mixture")
    mix_par <- names(mixture)
    mix_par_values <- mixture[[mix_par]]$values
    fits <- list()
    par_mix <- parameters
    ofvs <- c()
    for(i in seq(mix_par_values)) {
      par_mix[[mix_par]] <- mix_par_values[i]
      fits[[i]] <- mle_wrapper(
        minuslogl = ll_func,
        start = eta,
        method = method,
        optimizer = optimizer,
        control = control,
        skip_hessian = skip_hessian,
        data = list(
          data = data,
          sim_object = sim_object,
          parameters = par_mix,
          t_obs = t_obs,
          model = model,
          regimen = regimen,
          error = error,
          omega_full = omega_full_est / weight_prior,
          omega_inv = solve(omega_full_est / weight_prior),
          omega_eigen = sum(log(eigen(omega_full_est)$values)),
          nonfixed = nonfixed,
          transf = transf,
          weights = weights,
          weight_prior = weight_prior,
          sig = sig,
          as_eta = as_eta,
          censoring_idx = censoring_idx,
          censoring_label = censoring,
          iov_bins = iov_bins,
          calc_ofv = calc_ofv,
          covariates = covariates,
          steady_state_analytic = steady_state_analytic,
          include_omega = include_omega,
          include_error = include_error,
          verbose = verbose
        )
      )
      ofvs <- c(ofvs, fits[[i]]$log_likelihood)
    }
    like <- exp(ofvs)
    prior_prob <- c(mixture[[mix_par]]$probability, 1-mixture[[mix_par]]$probability)
    post_like <- prior_prob * like
    prob <- post_like / sum(post_like)
    mixture_group <- match(max(prob), prob)
    fit <- fits[[mixture_group]]
    parameters[[mix_par]] <- mix_par_values[mixture_group] # set population parameter value to most likely one
    mixture_obj <- list(
      parameter = mix_par,
      values = mix_par_values,
      mixture_group = mixture_group,
      selected = mix_par_values[mixture_group],
      probabilities = prob,
      prior_prob = prior_prob,
      likelihood = like
    )
  } else {
    output <- tryCatch({
      fit <- mle_wrapper(
        ll_func,
        start = eta,
        method = method,
        optimizer = optimizer,
        control = control,
        skip_hessian = skip_hessian,
        data = list(
          data = data,
          sim_object = sim_object,
          parameters = parameters,
          t_obs = t_obs,
          model = model,
          regimen = regimen,
          error = error,
          nonfixed = nonfixed,
          transf = transf,
          omega_full = omega_full_est / weight_prior,
          omega_inv = solve(omega_full_est / weight_prior),
          omega_eigen = sum(log(eigen(omega_full_est / weight_prior)$values)),
          weights = weights,
          weight_prior = weight_prior,
          sig = sig,
          as_eta = as_eta,
          censoring_idx = censoring_idx,
          censoring_label = censoring,
          iov_bins = iov_bins,
          calc_ofv = calc_ofv,
          covariates = covariates,
          steady_state_analytic = steady_state_analytic,
          include_omega = include_omega,
          include_error = include_error,
          verbose = verbose
        )
      )
    }, error = function(e) {
       return(e)
    })
    if("error" %in% class(output)) return(output)
  }
  cf <- fit$coef
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
      np <- get_np_estimates(
        parameter_grid = pars_grid,
        error = np_settings$error,
        model = model,
        regimen = regimen,
        data = data$y,
        t_obs = data$t,
        covariates = covariates,
        weights = weights
      )
      # take the estimates with highest probability as starting point for next grid
      tmp <- np$prob[order(-np$prob$like),][1,]
      par <- as.list(tmp[1:length(np$parameters)])
    }
    pars_grid <- create_grid_around_parameters(
      parameters = par,
      span = np_settings$grid_span,
      exponential = np_settings$grid_exponential,
      grid_size = np_settings$grid_size)
    np <- get_np_estimates(
      parameter_grid = pars_grid,
      error = np_settings$error,
      model = model,
      regimen = regimen,
      data = data$y,
      t_obs = data$t,
      covariates = covariates,
      weights = weights
    )
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
      if(!is.null(steady_state_analytic)) {
        A_init_ipred <- PKPDsim::calc_ss_analytic(
          f = steady_state_analytic$f,
          dose = regimen$dose_amts[1],
          interval = regimen$interval[1],
          model = model,
          parameters = obj$parameters,
          covariates = covariates,
          map = steady_state_analytic$map,
          n_transit_compartments = PKPDsim::ifelse0(steady_state_analytic$n_transit_compartments, FALSE),
          auc = PKPDsim::ifelse0(steady_state_analytic$auc, FALSE)
       )
       A_init_pred <- PKPDsim::calc_ss_analytic(
          f = steady_state_analytic$f,
          dose = regimen$dose_amts[1],
          interval = regimen$interval[1],
          model = model,
          parameters = parameters,
          covariates = covariates,
          map = steady_state_analytic$map,
          n_transit_compartments = PKPDsim::ifelse0(steady_state_analytic$n_transit_compartments, FALSE),
          auc = PKPDsim::ifelse0(steady_state_analytic$auc, FALSE)
       )
    } else {
      A_init_pred <- A_init
      A_init_ipred <- A_init
    }
    suppressMessages({
      sim_ipred <- PKPDsim::sim_ode(
        ode = model,
        parameters = par,
        mixture_group = mixture_group,
        covariates = covariates,
        n_ind = 1,
        int_step_size = int_step_size,
        regimen = regimen,
        t_obs = c(data_before_init$t, t_obs),
        obs_type = c(data_before_init$obs_type, data$obs_type),
        only_obs = TRUE,
        checks = FALSE,
        A_init = A_init_ipred,
        iov_bins = iov_bins,
        output_include = output_include,
        t_init = t_init,
        ...
      )
    })
    suppressMessages({
      sim_pred <- PKPDsim::sim_ode(
        ode = model,
        parameters = parameters,
        mixture_group = mixture_group,
        covariates = covariates,
        n_ind = 1,
        int_step_size = int_step_size,
        regimen = regimen,
        t_obs = c(data_before_init$t, t_obs),
        obs_type = c(data_before_init$obs_type, data$obs_type),
        only_obs = TRUE,
        checks = FALSE,
        iov_bins = iov_bins,
        A_init = A_init_pred,
        t_init = t_init,
        ...
      )
    })
    ipred <- sim_ipred$y
    pred <- sim_pred$y
    w_ipred <- sqrt(error$prop[data$obs_type]^2 * transf(ipred)^2 + error$add[data$obs_type]^2)
    w_pred <- sqrt(error$prop[data$obs_type]^2 * transf(pred)^2 + error$add[data$obs_type]^2)
    if(!(all(data$t == sim_ipred$t) && all(data$obs_type == sim_ipred$obs_type))) {
      warning("Mismatch in times and observation typese between input data and predictions. Be careful interpreting results from fit.")
    }
    y <- data$y
    prob <- list(par = c(mvtnorm::pmvnorm(cf, mean=rep(0, length(cf)),
                         sigma = omega_full[1:length(cf), 1:length(cf)])),
                 data = stats::pnorm(transf(y) - transf(ipred), mean = 0, sd = w_ipred))
    obj$res <- (transf(y) - transf(pred))
    obj$weights <- c(rep(0, length(data_before_init$t)), weights)
    obj$wres <- (obj$res / w_pred) * obj$weights
    obj$cwres <- obj$res / sqrt(abs(cov(transf(pred), transf(y_orig)))) * c(rep(0, nrow(data_before_init)), obj$weights)
    # Note: in NONMEM CWRES is on the population level, so can't really compare. NONMEM calls this CIWRES, it seems.
    obj$ires <- (transf(y_orig) - transf(ipred))
    obj$iwres <- (obj$ires / w_ipred)
    if(any(censoring_idx)) { # turn probabilities into "IWRES"-equivalent for censored data
      obj$iwres[censoring_idx] <- calc_res_from_prob(prob$data[censoring_idx])
      obj$ires[censoring_idx] <- obj$iwres[censoring_idx] * w_ipred[censoring_idx]
      ## if we would calculate the likelihood for the data population parameters given the data, 
      ## we could also calculate the the equivalents for cwres, wres, and res. However we 
      ## currently don't have a need to calculate. So setting to NA to avoid wrong interpretation.
      obj$cwres[censoring_idx] <- NA_real_
      obj$wres[censoring_idx] <- NA_real_
      obj$res[censoring_idx] <- NA_real_
    }
    obj$iwres_weighted <- obj$iwres * obj$weights
    obj$pred <- pred
    obj$ipred <- ipred
    obj$prob <- prob
    obj$dv <- y_orig
    obj$obs_type <- sim_ipred$obs_type
    if(output_include$covariates && !is.null(covariates)) {
      obj$covariates_time <- sim_ipred[!duplicated(sim_ipred$t), names(covariates)]
    }
    if(output_include$parameters) {
      obj$parameters_time <- sim_ipred[!duplicated(sim_ipred$t), names(parameters)]
    }
  }
  obj$vcov_full <- get_varcov_matrix(
    obj$fit$vcov, 
    fallback = omega_full
  )
  if(inherits(obj$vcov_full, "matrix")) {
    obj$vcov <- obj$vcov_full[t(!upper.tri(obj$vcov_full))]
  }
  obj$mahalanobis <- get_mahalanobis(y, ipred, w_ipred, ltbs)
  obj$prior <- list(
    parameters = parameters,
    omega = omega,
    fixed = fixed
  )
  class(obj) <- c(class(obj), "map_estimates")
  obj
}
