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
                      int_step_size = 0.25,
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
      if(nrow(data)>0) {
        data$evid <- 0
      } else {
        stop("No data available.")
      }
    }
  }
  colnames(data) <- tolower(colnames(data))
  if(!is.null(attr(model, "cpp")) && attr(model, "cpp")) {
    ll_func <- function(
      data,
      eta1 = 0, eta2 = 0, eta3 = 0, eta4 = 0, #eta5 = 0, eta6 = 0, eta7 = 0, eta8 = 0, eta9 = 0, eta10 = 0, eta11 = 0, eta12 = 0, # unfortunately seems no other way to do this...
      parameters,
      covariates,
      regimen = regimen,
      omega_full = omega_full,
      error = error,
      model,
      covs,
      min2LL = TRUE) {
        par <- parameters
        p <- as.list(match.call())
        for(i in 1:4) {
          par[[i]] <- par[[i]] * exp(p[[(paste0("eta", i))]])
        }
#        suppressMessages({
          sim <- sim_ode(ode = model,
                         parameters = par,
                         covariates = covariates,
                         n_ind = 1,
                         int_step_size = int_step_size,
                         # regimen = new_regimen(amt = dat[dat$evid == 1,]$amt, times = dat[dat$evid == 1,]$time, type = "infusion", t_inf = 2),
                         regimen = regimen,
                         t_obs = data[data$evid == 0,]$t) %>% dplyr::filter(comp == "obs")
          
#        })
        ipred <- sim[!duplicated(sim$t),]$y
        y <- data$y
        res_sd <- sqrt(error$prop^2*ipred^2 + error$add^2)
        ## need to adapt for different omega sizes!!
        ofv <-   c(dmvnorm(c(eta1, eta2, eta3, eta4), #, eta5, eta6, eta7, eta8, eta9, eta10, eta11, eta12), 
                           mean=rep(0, 4),
                           sigma=omega_full, 
                           log=TRUE),
                   dnorm(y - ipred, 0, res_sd, log=TRUE))
        if(verbose) { print(ofv) }
        if(min2LL) {
          return(-sum(ofv))
        } else {
          return(ofv)          
        }
      }
  } else {
    ll_func <- function(
      data,
      eta1 = 0, eta2 = 0, eta3 = 0, eta4 = 0, #, eta5 = 0, eta6 = 0, eta7 = 0, eta8 = 0, eta9 = 0, eta10 = 0, eta11 = 0, eta12 = 0, # unfortunately seems no other way to do this...
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
  for(i in 1:4) {
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
  if(lower_triangle_size(omega) < 4) {
    message("Omega matrix does not cover all parameters, will be expanded")
    add_length <- lower_triangle_vector_length(12) - length(omega) 
    omega <- c(omega, rep(0, add_length))
  }
  print(triangle_to_full(omega))
  print(eta)
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
#   res <- ll_func(data = data,
#                  eta1 = cf[1], 
#                  eta2 = cf[2], 
#                  parameters = parameters,
#                  covariates = covariates,
#                  regimen = regimen,
#                  omega_full = triangle_to_full(omega),
#                  error = error,
#                  model = model,
#                  covs = NULL,
#                  min2LL = FALSE)
  res <- NULL
  obj <- list(fit = fit, 
              parameters = par,
              residuals = res)
  class(obj) <- c(class(obj), "map_estimates")
  return(obj)
}
