% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/map_shrinkage_control.R
\name{map_shrinkage_control}
\alias{map_shrinkage_control}
\title{MAP fitting with shrinkage control}
\usage{
map_shrinkage_control(
  shrinkage_control = 0.05,
  parameters = NULL,
  fixed = NULL,
  data = NULL,
  ...
)
}
\arguments{
\item{shrinkage_control}{automatically control individual shrinkage. 
Suggested value e.g. `0.05` for allowing 5 percent shrinkage (averaged over 
all parameters).}

\item{parameters}{list of parameters}

\item{fixed}{fix a specific parameters, supplied as vector of strings}

\item{data}{data data.frame with columns `t` and `y` (and possibly evid)}

\item{...}{parameters passed on to `sim_ode()` function}
}
\description{
Limits the amount of individual shrinkage to a certain specified number
}
