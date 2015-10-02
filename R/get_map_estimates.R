#' @export
get_map_estimates <- function(
                      model = NULL,
                      data = NULL,
                      parameters = NULL,
                      covariates = NULL,
                      fixed = NULL,
                      omega = omega,
                      error = list(prop = 0.1, add = 0.1, exp = 0),
                      regimen = NULL,
                      int_step_size = 0.1,
                      method = "L-BFGS-B",
                      cols = list(x="t", y="y"),
                      verbose = FALSE) {
  if(is.null(model) || is.null(data) || is.null(parameters)) {
    stop("The 'model', 'data', and 'parameters' arguments are required.")
  }
  if(!("function" %in% class(model))) {
    stop("The 'model' argument requires a function, e.g. a model defined using the new_ode_model() function from the PKPDsim package.")
  }
  if(!all(unlist(cols) %in% names(data))) {
    stop("Expected column names were not found in data. Please use 'cols' argument to specify column names for independent and dependent variable.")
  }
  if("PKPDsim" %in% class(data)) {
    if("comp" %in% names(data)) {
      data <- data %>% dplyr::filter(comp == "obs")
      data <- data[!duplicated(data$t),]
      data$evid <- 0
    }
  }
  colnames(data) <- tolower(colnames(data))
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
      covs) {
        par <- parameters
        p <- as.list(match.call())
        for(i in seq(names(par))) {
          par[[i]] <- par[[i]] * exp(p[[(paste0("eta", i))]])
        }
        sig <- round(-log10(int_step_size))
        t_obs <- round(data[data$evid == 0,]$t, sig)
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
        ofv <-   c(dmvnorm(c(eta1, eta2), mean=c(0, 0), sigma=omega_full, log=TRUE),
                   dnorm(y - ipred, 0, res_sd, log=TRUE))
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
                   dnorm(y - ipred, 0, res_sd, log=TRUE))
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
  fit <- mle2(ll_func,
              start = eta,
              method = method,
              data = list(data = data,
                          parameters = parameters,
                          covariates = covariates,
                          regimen = regimen,
                          model = model,
                          omega_full = triangle_to_full(omega),
                          error = error,
                          covs = NULL),
              fixed = fix)
  cf <- coef(fit)
  par <- parameters
  for(i in seq(names(par))) {
    par[[i]] <- as.numeric(as.numeric(par[[i]]) * exp(as.numeric(cf[i])))
  }
  obj <- list(fit = fit, parameters = par)
  class(obj) <- c(class(obj), "map_estimates")
  return(obj)
}
