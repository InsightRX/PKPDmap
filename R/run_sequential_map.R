#' Perform a sequential MAP Bayesian fit.
#'
#' This will allow "tracking" parameters over time, i.e. a form of between-occasion variability
#'
#' @param model PKPDsim model
#' @param regimen PKPDsim regimen
#' @param parameters population values for parameters for `model``
#' @param covariates covariates for `model`, if present
#' @param data data.frame with `t` and `y` columns
#' @param breaks_by_data set breaks automatically (`T`/`F``), 1 data point per section (end of section).
#' @param breaks_by_time set breaks automatically by time (`T`/`F``), also specify `break_interval`
#' @param break_interval when `breaks_by_time` is `TRUE`, the interval between breaks
#' @param breaks vector of breaks for sections
#' @param A_init vector of initial state
#' @param weight_prior weight of population priors specified through `omega`
#' @param weight_unfocus weight for the points in the sequential fit that are not in the focused section (but later).
#' @param ... additional arguments passed to `get_map_estimates()` and `sim_ode()`.
#' @param verbose verbose output?
#' @export
run_sequential_map <- function(
  model = NULL,
  regimen = NULL,
  parameters = list(),
  covariates = NULL,
  data = NULL,
  breaks_by_data = TRUE,
  breaks_by_time = FALSE,
  break_interval = 24,
  breaks = NULL,
  weights_prior = NULL,
  A_init = NULL,
  weight_prior = 1,
  weight_unfocus = .1,
  verbose = FALSE,
  ...
) {
  weights_prior <- 0.5
  fits <- list()
  dats <- list()
  if(is.null(data)) {
    stop("No data supplied to fit!")
  }
  if(is.null(A_init)) {
    A_init <- rep(0, attr(model, "size"))
  }
  conc <- c()
  par_table <- c()
  if(is.null(breaks)) {
    if(breaks_by_data) {
      breaks <- data$t
    }
    if(breaks_by_time) {
      breaks <- unique(c(seq(from = 0, to = max(data$t), by = break_interval), max(data$t)))
    }
  }
  breaks <- unique(c(0, breaks))

  ## test if all sections have observations, which is req'd for MAP
  tmp <- cut(data$t, breaks)
  if(!all(levels(tmp) %in% unique(tmp[!is.na(tmp)]))) {
    filt <- levels(tmp) %in% unique(tmp[!is.na(tmp)])
    last_break <- breaks[1]
    for(i in 1:length(filt)) {
      if(!filt[i]) {
        breaks[i+1] <- last_break
      }
      last_break <- breaks[i+1]
    }
    breaks <- unique(breaks)
    cat(paste0("Note: all sections need to contain at least one datapoint, merging empty sections to: ", paste(breaks, collapse=" "), ".\n\n"))
    ## merge sections
  }
  if(verbose) {
    pb <- txtProgressBar(min = 1, max = length(breaks))
  }
  t_max_calc <- regimen$dose_times + regimen$interval
  for(i in 2:length(breaks)) {
    if(verbose) { setTxtProgressBar(pb, i) }
    reg_tmp <- regimen
    obs_tmp <- data
    t_last <- 0
    if(i > 1) { # update data
      t_last <- breaks[i-1]
      reg_tmp$dose_times <- reg_tmp$dose_times - t_last
      reg_tmp$dose_times <- reg_tmp$dose_times[reg_tmp$dose_times >= 0]
      reg_tmp$dose_amts <- tail(reg_tmp$dose_amts, length(reg_tmp$dose_times))
      reg_tmp$t_inf <- tail(reg_tmp$t_inf, length(reg_tmp$dose_times))
      reg_tmp$type <- tail(reg_tmp$type, length(reg_tmp$dose_times))
      obs_tmp <- data[data$t > t_last,]
      obs_tmp$t <- obs_tmp$t - t_last
      # covariates too
    }
    t1 <- 0
    t2 <- breaks[i] - t_last
    weights <- rep(weight_unfocus, length(obs_tmp$t))
    weights[obs_tmp$t < t2] <- 1

    ## Fit to data from i to end
    fits[[i]] <- PKPDmap::get_map_estimates(
      model = model,
      data = obs_tmp,
      A_init = A_init,
      regimen = reg_tmp,
      omega = omega,
      covariates = covariates,
      weights = weights,
      weight_prior = weight_prior,
      parameters = parameters,
      ruv = ruv,
      verbose = verbose)
    ## simulate out data from i-1 to i
    tmp <- PKPDsim::sim_ode(ode = model,
                   regimen = reg_tmp,
                   covariates = covariates,
                   parameters = fits[[i]]$parameters,
                   A_init = A_init,
                   t_obs=c(seq(from = t1, to = floor(t2), by = 0.1), t2),
                   verbose = verbose)
    conc <- rbind(conc,
                  tmp %>% filter(comp == "obs") %>% mutate(t = t + t_last))
    par_table <- rbind(par_table,
                       cbind(t1 = t1 + t_last,
                             t2 = t2 + t_last, as.data.frame(fits[[i]]$parameters)))
    tmp_filt <- tmp %>%
      dplyr::filter(t == max(t) & !duplicated(t) & comp != "obs")
    A_init <- tmp_filt$y
  }
  class(conc) <- c("PKPDsim_data", "data.frame")
  return(list(
    fit = fits,
    obs = conc,
    parameters = par_table
  ))
}
