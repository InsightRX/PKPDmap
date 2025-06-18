#' Parse input data for fit
#' 
#' @inheritParams get_map_estimates
#' 
parse_input_data <- function(data) {
  if(inherits(data, "PKPDsim")) {
    # when fitting directly from PKPDsim-generated data:
    if("comp" %in% names(data)) {
      data <- data[data$comp == "obs",]
      data$evid <- 0
    }
  }
  colnames(data) <- tolower(colnames(data))
  if("evid" %in% colnames(data)) {
    data <- data[data$evid == 0,]
  }
  data[order(data$t, data$obs_type),]
}