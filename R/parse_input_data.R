#' Parse input data for fit
#' 
#' @inheritParams get_map_estimates
#' 
parse_input_data <- function(data, obs_type_label = NULL) {
  colnames(data) <- tolower(colnames(data))
  if(inherits(data, "PKPDsim")) {
    # when fitting directly from PKPDsim-generated data:
    if("comp" %in% names(data)) {
      data <- data[data$comp == "obs",]
      data$evid <- 0
    }
  }
  ## Parse data to include label for obs_type
  if(is.null(obs_type_label)) {
    data$obs_type <- 1
  } else {
    data$obs_type <- data[[obs_type_label]]
  }
  if("evid" %in% colnames(data)) {
    data <- data[data$evid == 0,]
  }
  data[order(data$t, data$obs_type),]
}