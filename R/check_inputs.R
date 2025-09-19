#' Some basic input checking
#' 
#' @inheritParams get_map_estimates
#' 
check_inputs <- function(model, data, parameters, omega, regimen, censoring, type) {
  if(tolower(type) %in% c("map", "pls")) {
    if(is.null(model) || is.null(data) || is.null(parameters) || is.null(omega) || is.null(regimen)) {
      stop("The 'model', 'data', 'omega', 'regimen', and 'parameters' arguments are required.")
    }
  }
  if(!is.null(censoring) && !inherits(censoring, "character")) {
    stop("Censoring argument requires label specifying column in dataset with censoring info.")
  }
  if(!("function" %in% class(model))) {
    stop("The 'model' argument requires a function, e.g. a model defined using the new_ode_model() function from the PKPDsim package.")
  }
  if(!is.null(parameters)) {
    defined_parameters <- names(attr(model, "parameters"))
    if(! all(defined_parameters %in% names(parameters))) {
      stop("One or more required parameters for the model have not been specified.")
    }
    if(any(names(parameters) %in% defined_parameters)) {
      warning("One or more of the provided `parameters` are not supported by the model and will be ignored. Passing unknown parameters may affect IIV and IOV structure and result in erroneous output.")
    }
  }
}