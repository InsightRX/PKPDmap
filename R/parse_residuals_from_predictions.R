#' Parse simulated data and extract residuals 
#' 
#' @inheritParams calc_residuals
#' @param sim_ipred data frame with simulated individual predictions
#' @param sim_pred data frame with simulated population predictions
#' 
parse_residuals_from_predictions <- function(
  obj,
  sim_ipred,
  sim_pred,
  data,
  omega_full,
  transf,
  error,
  censoring,
  censoring_idx,
  data_before_init,
  weights
) {
  ipred <- sim_ipred$y
  pred <- sim_pred$y
  w_ipred <- sqrt(error$prop[data$obs_type]^2 * transf(ipred)^2 + error$add[data$obs_type]^2)
  w_pred <- sqrt(error$prop[data$obs_type]^2 * transf(pred)^2 + error$add[data$obs_type]^2)
  if(!(all(data$t == sim_ipred$t) && all(data$obs_type == sim_ipred$obs_type))) {
    warning("Mismatch in times and observation typese between input data and predictions. Be careful interpreting results from fit.")
  }
  y <- data$y
  cf <- obj$fit$coef
  prob <- list(
    par = c(
      mvtnorm::pmvnorm(
        cf,
        mean=rep(0, length(cf)),
        sigma = omega_full[1:length(cf), 1:length(cf)]
      )
    ),
    data = stats::pnorm(transf(y) - transf(ipred), mean = 0, sd = w_ipred)
  )
  obj$res <- (transf(y) - transf(pred))
  obj$weights <- c(rep(0, length(data_before_init$t)), weights)
  obj$wres <- (obj$res / w_pred) * obj$weights
  obj$cwres <- obj$res / sqrt(abs(cov(transf(pred), transf(y)))) * c(rep(0, nrow(data_before_init)), obj$weights)
  # Note: in NONMEM CWRES is on the population level, so can't really compare. NONMEM calls this CIWRES, it seems.
  obj$ires <- (transf(y) - transf(ipred))
  obj$iwres <- (obj$ires / w_ipred)
  obj$w_ipred  <- w_ipred
  if(is.null(censoring)) {
    obj$censoring <- rep(0, length(y))
  } else {
    obj$censoring <- data[[censoring]]
    if(any(censoring_idx)) {
      obj$iwres[censoring_idx] <- calc_res_from_prob(prob$data[censoring_idx])
      obj$ires[censoring_idx] <- obj$iwres[censoring_idx] * w_ipred[censoring_idx]
      ## if we would calculate the likelihood for the data population parameters given the data, 
      ## we could also calculate the the equivalents for cwres, wres, and res. However we 
      ## currently don't have a need to calculate. So setting to NA to avoid wrong interpretation.
      obj$cwres[censoring_idx] <- NA_real_
      obj$wres[censoring_idx] <- NA_real_
      obj$res[censoring_idx] <- NA_real_
    }
  }
  obj$iwres_weighted <- obj$iwres * obj$weights
  obj$pred <- pred
  obj$ipred <- ipred
  obj$prob <- prob
  obj$dv <- y
  obj$obs_type <- sim_ipred$obs_type
  obj
}