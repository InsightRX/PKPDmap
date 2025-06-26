#' Check if censoring functionality needs to be used, 
#' and create index for censored data
#' 
#' @inheritParams get_map_estimates
#'  
parse_censoring <- function(
  censoring, 
  data,
  verbose = FALSE
) {
  censoring_idx <- NULL
  if(!is.null(censoring)) {
    if(any(data[[tolower(censoring)]] != 0)) {
      censoring_idx <- data[[tolower(censoring)]] != 0
      if(verbose) message("One or more values in data are censored, including censoring in likelihood.")
    } else {
      if(verbose) message("Warning: censoring specified, but no censored values in data.")
    }
  }
  censoring_idx
}