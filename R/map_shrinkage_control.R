#' MAP fitting with shrinkage control
#' 
#' Limits the amount of individual shrinkage to a certain specified number
#'
#' @inheritParams get_map_estimates
#' 
#' @param shrinkage_control automatically control individual shrinkage. 
#' Suggested value e.g. `0.05` for allowing 5 percent shrinkage (averaged over 
#' all parameters).
#' 
#' @export
#' 
map_shrinkage_control <- function(
  shrinkage_control = 0.05,
  parameters = NULL,
  fixed = NULL,
  data = NULL,
  ...
) {
  weights <- seq(0.1, 1, 0.1)
  fits <- c()
  fit_list <- list()
  for(i in seq(weights)) {
    fit_list[[i]] <- get_map_estimates(
      weight_prior = weights[i],
      parameters = parameters,
      fixed = fixed,
      data = data,
      ...)
    fits <- rbind(fits, c(weights[i], unlist(fit_list[[i]]$parameters)))
  }
  nonfixed <- names(parameters)[is.na(match(names(parameters), fixed))]
  fits <- data.frame(fits)
  names(fits) <- c("weight", nonfixed)
  for(i in seq(nonfixed)) {
    p <- nonfixed[i]
    regr <- lm(data = fits, formula = get(p) ~ poly(weight,2))
    zero_shrink <- as.numeric(predict.lm(regr, data.frame(weight=0)))
    fits[[paste(p, "_shr")]] <- abs((fits[[p]] - zero_shrink) / (zero_shrink - parameters[[p]]))
  }
  # pick fit closest to desired average shrinkage
  shr_cols <- grep("_shr", names(fits))
  closest <- order(abs(apply((fits[,shr_cols]) - shrinkage_control, 1, mean)))[1]
  fit <- fit_list[[closest]]
  fit$shrinkage_control <- list(
    weight_prior = weights[closest],
    fit_list = fit_list
  )
  return(fit)
}
