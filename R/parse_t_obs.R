#' Parse observation times
#' 
#' @inheritParams get_map_estimates
#' 
parse_t_obs <- function(data) {
  t_obs <- data$t
  if(any(duplicated(paste(t_obs, data$obs_type, sep = "_")))) {
    message("Duplicate times were detected in data. Estimation will proceed but please check that data is correct. For putting more weight on certain measurements, please use the `weights` argument.")
  }
  t_obs
}
