#' Print MAP estimates
#' 
#' @param x output object from get_map_estimates
#' @param ... etc
#' @export
print.map_estimates <- function(x, ...) {
  eta <- x$fit$coef
  om <- sqrt(diag(PKPDsim::triangle_to_full(x$prior$omega)))
  res_table <- data.frame(
    individual = unlist(x$parameters), 
    population = unlist(x$prior$parameters)
  )
  res_table$eta[! names(x$parameters) %in% x$prior$fixed] <- eta / om[1:length(eta)]
  res_table$relative[! names(x$parameters) %in% x$prior$fixed] <- plot_eta(eta / om[1:length(eta)])
  print(res_table)
  gof <- data.frame(
    dv = x$dv,
    ipred = x$ipred,
    pred = x$pred,
    res = x$res,
    iwres = x$iwres
  )
  if(length(unique(x$obs_type)) > 1) {
    gof$obs_type <- x$obs_type
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
  if(length(x) > 1) {
    return(unlist(lapply(x, "plot_eta")))
  }
  dummy <- "-----"
  if(x == 0) {
    return(paste0(dummy, "o", dummy))
  }
  position <- min(round(abs(x)) * 2, 5) # the factor 2 is arbitrary to have a nice scale.
  scale <- paste0(
    paste0(rep("-", position), collapse=""), 
    "o", 
    paste0(rep("-", max(4-position, 0)), collapse=""), 
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
