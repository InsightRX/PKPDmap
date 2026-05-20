#' Calculate CWRES (Conditional Weighted Residuals)
#'
#' Computes CWRES following the FOCE approximation method described in
#' Hooker et al. (2007). The Jacobian of model predictions with respect to
#' random effects (etas) is computed via central finite differences.
#'
#' CWRES = L^{-1} * (y - f(eta_hat) + F * eta_hat)
#'
#' where:
#' - f(eta_hat) = individual predictions (IPRED) in the transformed space
#' - F = Jacobian df/deta evaluated at eta_hat
#' - L * L^T = V = F * Omega * F^T + Sigma (FOCE variance)
#' - Sigma = diagonal residual error variance matrix
#'
#' @references Hooker AC, Staatz CE, Karlsson MO. Conditional weighted
#'   residuals (CWRES): a model diagnostic for the FOCE method.
#'   Pharm Res. 2007;24(12):2187-2197.
#'
#' @param eta_hat numeric vector of estimated etas (MAP estimates)
#' @param ipred_raw numeric vector of raw (untransformed) individual
#'   predictions at eta_hat
#' @param y numeric vector of observed values
#' @param omega_full full omega covariance matrix
#' @param error list with `prop` and `add` residual error components
#' @param obs_type integer vector of observation types
#' @param transf transformation function (identity or log for LTBS)
#' @param model PKPDsim model object
#' @param parameters_population list of population parameters
#' @param nonfixed character vector of non-fixed parameter names
#' @param as_eta character vector of parameters estimated directly as eta
#' @param covariates list of PKPDsim covariates
#' @param regimen PKPDsim regimen object
#' @param lagtime lagtime specification
#' @param t_obs numeric vector of observation times
#' @param obs_type_sim observation type vector for simulation
#' @param int_step_size integrator step size
#' @param iov_bins IOV bin specification
#' @param A_init initial state vector
#' @param t_init initialization time
#' @param steady_state_analytic steady state settings (or NULL)
#' @param weight_prior_var prior weight variance scaling factor (default 1).
#'   Used to scale omega in the FOCE Hessian computation to match the
#'   estimation objective. CWRES uses the unscaled omega (model diagnostic).
#' @param delta perturbation size for finite differences
#'
#' @return list with components:
#'   \item{cwres}{numeric vector of CWRES values}
#'   \item{vcov}{FOCE variance-covariance matrix of eta estimates
#'     (n_eta x n_eta), or NULL if computation failed. Derived from the
#'     same Jacobian F used for CWRES: vcov = (Omega^-1 + F' Sigma^-1 F)^-1}
#'
calc_cwres <- function(
    eta_hat,
    ipred_raw,
    y,
    omega_full,
    error,
    obs_type,
    transf,
    model,
    parameters_population,
    nonfixed,
    as_eta,
    covariates,
    regimen,
    lagtime,
    t_obs,
    obs_type_sim,
    int_step_size = 0.1,
    iov_bins = NULL,
    A_init = NULL,
    t_init = 0,
    steady_state_analytic = NULL,
    weight_prior_var = 1,
    delta = 1e-4
) {
  n_obs <- length(y)
  n_eta <- length(eta_hat)

  if (n_obs == 0 || n_eta == 0) {
    return(list(cwres = numeric(0), vcov = NULL))
  }

  # Transformed individual predictions at eta_hat
  ipred_transf <- transf(ipred_raw)

  # Compute Jacobian F = d(transf(f))/d(eta) via central finite differences
  F_matrix <- matrix(0, nrow = n_obs, ncol = n_eta)
  for (j in seq_len(n_eta)) {
    eta_plus <- eta_hat
    eta_minus <- eta_hat
    eta_plus[j] <- eta_hat[j] + delta
    eta_minus[j] <- eta_hat[j] - delta

    pred_plus <- simulate_with_etas(
      eta_plus, parameters_population, nonfixed, as_eta,
      model, covariates, regimen, lagtime, t_obs, obs_type_sim,
      int_step_size, iov_bins, A_init, t_init,
      steady_state_analytic
    )
    pred_minus <- simulate_with_etas(
      eta_minus, parameters_population, nonfixed, as_eta,
      model, covariates, regimen, lagtime, t_obs, obs_type_sim,
      int_step_size, iov_bins, A_init, t_init,
      steady_state_analytic
    )

    F_matrix[, j] <- (transf(pred_plus) - transf(pred_minus)) / (2 * delta)
  }

  # Residual error variance diagonal (in transformed space)
  sigma_diag <- error$prop[obs_type]^2 * ipred_transf^2 +
    error$add[obs_type]^2

  if (any(!is.finite(sigma_diag) | sigma_diag <= 0)) {
    bad <- which(!is.finite(sigma_diag) | sigma_diag <= 0)
    warning(
      "CWRES/vcov computation failed: residual error variance is zero, ",
      "negative, or non-finite for observation(s) ",
      paste(bad, collapse = ", "),
      " (obs_type: ", paste(unique(obs_type[bad]), collapse = ", "), ")."
    )
    return(list(cwres = rep(NA_real_, n_obs), vcov = NULL))
  }

  # Omega submatrix for estimated parameters
  omega_est <- omega_full[seq_len(n_eta), seq_len(n_eta), drop = FALSE]

  # FOCE approximate marginal variance: V = F * Omega * F^T + Sigma
  V <- F_matrix %*% omega_est %*% t(F_matrix) +
    diag(sigma_diag, nrow = n_obs)

  # Cholesky decomposition (lower triangular)
  L <- tryCatch(
    t(chol(V)),
    error = function(e) NULL
  )
  cwres <- if (is.null(L)) {
    warning("CWRES computation failed: variance matrix is not positive definite.")
    rep(NA_real_, n_obs)
  } else {
    # CWRES = L^{-1} * (y_transf - ipred_transf + F * eta_hat)
    # Derivation:
    #   E_FOCE(y) = f(eta_hat) - F * eta_hat  (since E[eta] = 0)
    #   y - E_FOCE(y) = y - f(eta_hat) + F * eta_hat
    y_transf <- transf(y)
    as.numeric(
      solve(L, y_transf - ipred_transf + F_matrix %*% eta_hat)
    )
  }

  # FOCE Hessian: H = (Omega/w)^{-1} + F' * Sigma^{-1} * F
  # where w = weight_prior_var, matching the scaled omega used during MAP
  # estimation. vcov of etas = H^{-1}.
  # This reuses the Jacobian F already computed for CWRES, so no extra
  # simulations are needed (replaces the numDeriv::hessian computation).
  # With weight_prior_var <= 0 the prior is effectively flat (LS fit), so
  # the FOCE vcov is not meaningful; skip and let the caller fall back.
  vcov <- NULL
  if (weight_prior_var > 0) {
    omega_scaled <- omega_est / weight_prior_var
    omega_scaled_inv <- tryCatch(solve(omega_scaled), error = function(e) NULL)
    if (!is.null(omega_scaled_inv)) {
      H_foce <- omega_scaled_inv +
        t(F_matrix) %*% diag(1 / sigma_diag, nrow = n_obs) %*% F_matrix
      vcov <- tryCatch(solve(H_foce), error = function(e) {
        warning("FOCE variance-covariance computation failed.")
        NULL
      })
    }
  }

  list(cwres = cwres, vcov = vcov)
}


#' Simulate predictions with a given eta vector
#'
#' Helper for computing the Jacobian in CWRES calculation via finite
#' differences. Computes individual parameters from etas and runs a
#' PKPDsim simulation.
#'
#' @keywords internal
simulate_with_etas <- function(
    eta,
    parameters_population,
    nonfixed,
    as_eta,
    model,
    covariates,
    regimen,
    lagtime,
    t_obs,
    obs_type,
    int_step_size,
    iov_bins,
    A_init,
    t_init,
    steady_state_analytic
) {
  # Compute individual parameters from etas
  par <- parameters_population
  for (i in seq_along(nonfixed)) {
    key <- nonfixed[i]
    if (key %in% as_eta) {
      par[[key]] <- eta[i]
    } else {
      par[[key]] <- par[[key]] * exp(eta[i])
    }
  }

  # Recompute steady-state A_init if needed
  a_init <- A_init
  if (!is.null(steady_state_analytic)) {
    a_init <- PKPDsim::calc_ss_analytic(
      f = steady_state_analytic$f,
      dose = regimen$dose_amts[1],
      interval = regimen$interval[1],
      model = model,
      parameters = par,
      covariates = covariates,
      map = steady_state_analytic$map,
      n_transit_compartments = PKPDsim::ifelse0(
        steady_state_analytic$n_transit_compartments, FALSE
      ),
      auc = PKPDsim::ifelse0(steady_state_analytic$auc, FALSE)
    )
  }

  suppressMessages({
    sim <- PKPDsim::sim_ode(
      ode = model,
      parameters = par,
      mixture_group = NULL,
      covariates = covariates,
      n_ind = 1,
      int_step_size = int_step_size,
      regimen = regimen,
      t_obs = t_obs,
      obs_type = obs_type,
      only_obs = TRUE,
      checks = FALSE,
      A_init = a_init,
      iov_bins = iov_bins,
      t_init = t_init,
      lagtime = lagtime
    )
  })

  sim$y
}
