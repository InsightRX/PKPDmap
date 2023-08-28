#' Calculate residuals from probabilities for datapoints
#' 
#' Calculates the residual associated with a certain probability
#' It's basically the inverse of what we do in the likelihood function. There
#' we calculate the probability of a certain point, by calculating
#' 
#'   p = dnorm(obs - ipred, mean, sd) / dnorm(0, mean, sd).
#' 
#' Here, we want to do the reverse, i.e. we have a probability but want to get
#' the residual. Reshuffling the above equation and taking the inverse of dnorm
#' gives us that. This is useful for censored observations, where we don't 
#' actually have the "observed" value, but we can still calculate a residual-
#' equivalent using this function.
#'
#' Note: we actually get back the abs(residual), since from only the probability
#' we cannot infer the direction of the residual.
#' 
#' @param p probability, e.g. from fit object likelihood info
calc_res_from_prob <- function(p) {
  p <- pmin(p, 0.9999) # avoid Inf at 1
  p <- pmax(p, 0.0001) # avoid Inf at 0
  dnorminv(
    p * dnorm(0, mean = 0, sd = 1)
  )
}

#' Inverse of the dnorm function in R. 
#' 
#' This is not the same as qnorm, which is the inverse of pnorm (the CDF).
#' 
#' Taken from: https://stackoverflow.com/questions/19589191/the-reverse-inverse-of-the-normal-distribution-function-in-r
#' 
#' @param y density
dnorminv <- function(y) {
  sqrt(-2*log(sqrt(2*pi)*y))
}