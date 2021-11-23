#' Print MAP estimates
#' 
#' @param x output object from get_map_estimates
#' @param ... etc
#' @export
print.map_estimates <- function(x, ...) {
  print(data.frame(
    individual = unlist(x$parameters), 
    population = unlist(x$prior$parameters),
    eta = x$fit$coef,
    location = unlist(lapply(x$fit$coef, plot_eta))
  ))
  gof <- data.frame(
    dv = fit$dv,
    ipred = fit$ipred,
    pred = fit$pred,
    res = fit$res,
    iwres = fit$iwres
  )
  if(length(unique(fit$obs_type)) > 1) {
    gof$obs_type <- fit$obs_type
  }
  cat("\n")
  print(gof)
  cat(paste0("\nOFV:\t", signif(x$fit$value, 6)))
}

#' Visual representation of eta estimate
#' 
#' @param x estimated individual coefficient for parameter (=eta in NONMEM) 
#' @export
plot_eta <- function(x) {
  dummy <- "-----"
  position <- min(round(abs(x) * 4), 5)
  scale <- paste0(
    paste0(rep("-", position), collapse=""), 
    "o", 
    paste0(rep("-", 4-position), collapse=""), 
    collapse = "")
  out <- paste0(
    dummy, 
    "|",
    scale,
    collapse = ""
  )
  if(x < 0) { # reverse
    splits <- strsplit(out, "")[[1]]
    reversed <- rev(splits)
    out <- paste(reversed, collapse = "")
  }
  out
}