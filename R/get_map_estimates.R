#' Get MAP estimates
#' 
#' @param model model
#' @param data data
#' @param parameters parameters
#' @param covariates covariates
#' @param fixed fix a specific parameters, supplied as vector of strings
#' @param omega between subject variability, supplied as vector specifiying the lower triangle of the covariance matrix of random effects
#' @param error residual error, specified as list with arguments `add` and/or `prop` specifying the additive and proportional parts 
#' @param regimen regimen
#' @param int_step_size integrator step size passed to PKPDsim
#' @param method optimization method, default L-BFGS-B
#' @param cols column names
#' @param verbose show more output
#' @export
get_map_estimates <- function(
                      model = NULL,
                      data = NULL,
                      parameters = NULL,
                      covariates = NULL,
                      fixed = NULL,
                      omega = NULL,
                      error = list(prop = 0.1, add = 0.1, exp = 0),
                      regimen = NULL,
                      int_step_size = 0.1,
                      method = "L-BFGS-B",
                      cols = list(x="t", y="y"),
                      verbose = FALSE) {
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
  t_obs <- data$t
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
                         only_obs = TRUE)
        })
        ipred <- sim[!duplicated(sim$t),]$y
        y <- data$y
        res_sd <- sqrt(error$prop^2*ipred^2 + error$add^2)
        ## need to adapt for different omega sizes!!
        ofv <-   c(mvtnorm::dmvnorm(c(eta1, eta2), mean=c(0, 0), sigma=omega_full, log=TRUE),
                   dnorm(y - ipred, mean = 0, sd = res_sd, log=TRUE))
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
        ## need to adapt for different omega sizes!!
        ofv <-   c(dmvnorm(c(eta1, eta2), mean=c(0, 0), sigma=omega_full, log=TRUE),
                   dnorm(y - ipred, 0, sd = res_sd, log=TRUE))
        if(verbose) { print(ofv) }
        return(-sum(ofv))
      }
  }
  eta <- list()
  for(i in seq(names(parameters))) {
    eta[[paste0("eta", i)]] <- 0
  }
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
  fit <- bbmle::mle2(ll_func,
              start = eta,
              method = method,
              data = list(data = data,
                          parameters = parameters,
                          covariates = covariates,
                          regimen = regimen,
                          model = model,
                          omega_full = triangle_to_full(omega),
                          error = error,
                          t_obs = t_obs,
                          sig = sig,
                          covs = NULL),
              fixed = fix)
  cf <- bbmle::coef(fit)
  par <- parameters
  for(i in seq(names(par))) {
    par[[i]] <- as.numeric(as.numeric(par[[i]]) * exp(as.numeric(cf[i])))
  }
  obj <- list(fit = fit, parameters = par)
  class(obj) <- c(class(obj), "map_estimates")
  return(obj)
}
