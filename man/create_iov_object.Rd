% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_iov_object.R
\name{create_iov_object}
\alias{create_iov_object}
\title{Calculate the specs for including IOV, e.g. which elements should be fixed
and the right size of the new omega matrix.
Currently only working for a single parameter!}
\usage{
create_iov_object(
  cv = list(CL = 0.1),
  omega = c(0.1),
  bins = c(0, 24, 48, 9999),
  parameters = list(CL = 5),
  tdm_init_level = NULL,
  estimate_init_level = FALSE,
  init_level_weight = 0.5,
  ruv = NULL,
  fixed = NULL,
  n = NULL,
  verbose = TRUE
)
}
\arguments{
\item{cv}{list of CVs of parameters}

\item{omega}{omega matrix (lower triangle)}

\item{bins}{vector of bins for IOV}

\item{parameters}{named list of parameter values}

\item{tdm_init_level}{pre-first dose TDM level}

\item{estimate_init_level}{estimate the pre-dose. If FALSE, will use any TDM levels before first dose as the deterministic level for compartment initiation.}

\item{init_level_weight}{weight in MAP fit for init level}

\item{ruv}{residual variability model. Required when estimate_init_level is `TRUE` to estimate error on init level.}

\item{fixed}{vector of fixed parameters}

\item{n}{number of IOV elements, will be determined automatically from data, `n` will override if not NULL.}

\item{verbose}{verbosity (`TRUE` or `FALSE`)}
}
\description{
Calculate the specs for including IOV, e.g. which elements should be fixed
and the right size of the new omega matrix.
Currently only working for a single parameter!
}
