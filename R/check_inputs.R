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
  check_parameters_matching(model, parameters)
}

check_parameters_matching <- function(model, parameters) {
  defined_parameters <- attr(model, "parameters")
  if(is.null(defined_parameters)) {
    warning("Parameter information for model missing, cannot perform parameter consistency check. Please check PKPDsim model definition.")
  } else {
    if(! all(defined_parameters %in% names(parameters))) {
      missing_pars <- defined_parameters[! defined_parameters %in% names(parameters)]
      stop(
        paste0(
          "One or more required parameters for the model have not been specified. Missing: ",
          paste0(missing_pars, collapse = ", ")
        )
      )
    }
    if(any(!(names(parameters) %in% defined_parameters))) {
      ignored_pars <- names(parameters)[! names(parameters) %in% defined_parameters]
      warning(
        paste0(
          "Some supplied `parameters` are not supported by the model and will be ignored: ",
          paste0(ignored_pars, collapse = ", "),
          ". Passing unknown parameters may affect IIV and IOV structure and result in erroneous output."
        )
      )
    }
  }
}
