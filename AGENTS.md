# AGENTS.md

This file provides guidance to AI coding agents when working with code in this
repository. It is intentionally harness-agnostic: any agent tooling that reads
`AGENTS.md` should be able to use it as-is.

## Package Overview

PKPDmap is an R package implementing Maximum A Posteriori (MAP) Bayesian
estimation for pharmacokinetic/pharmacodynamic (PK/PD) data. It depends heavily
on [PKPDsim](https://github.com/InsightRX/PKPDsim) (also an InsightRX package)
for ODE-based PK/PD simulation. PKPDsim is not on CRAN and is pulled from
GitHub via the `Remotes:` field in `DESCRIPTION`.

## Common Commands

All development uses standard R/devtools workflows (no Makefile):

```r
# Load package for interactive development
devtools::load_all()

# Run all tests
devtools::test()

# Run a single test file (matches on the part after "test-")
devtools::test(filter = "get_map_estimates")

# Regenerate documentation from roxygen2 comments
devtools::document()

# Full R CMD check (same flags as CI)
devtools::check(args = c("--no-manual", "--as-cran"))

# Install PKPDsim dependency from GitHub (public repo; no token needed, but
# `GITHUB_PAT` can be set to avoid API rate limits)
remotes::install_github("InsightRX/PKPDsim")
```

## Architecture

### Estimation Types

The main entry point is `get_map_estimates()` (`R/get_map_estimates.R`). Two
separate arguments are easy to confuse:

- **`type`** selects the *estimation approach*. `check_inputs()` and
  `parse_weight_prior()` recognize:
  - `"map"` (default) — standard MAP Bayesian / empirical Bayes. Uses
    `calc_ofv_map()`, which adds the prior penalty term.
  - `"pls"` — penalized least squares: MAP with very flat priors. Implemented
    by forcing `weight_prior_var` to `0.001` (empirically chosen).
  - `"ls"` — least squares: switches to `calc_ofv_ls()` and overrides the
    residual error model to `list(prop = 0, add = 1)`.
- **`method`** is passed straight through to the optimizer (default `"BFGS"`).
  It is `optim()`'s method, *not* the estimation type. BFGS is gradient-based
  with finite-difference gradients; for models containing step functions (e.g.
  a lagtime) prefer a non-gradient method such as Nelder-Mead.

Note: `README.md` documents these as `map` / `map_flat_prior` / `ls`. The code
actually accepts `map` / `pls` / `ls` — treat the source as authoritative.

Prior strength is controlled independently by **`weight_prior`**, which is
given on the SD scale and squared into `weight_prior_var`. Throughout the
pipeline the prior is applied as `omega$est / weight_prior_var`, so a larger
`weight_prior` means a *tighter* prior. Do **not** use `weight_prior = 0` to
drop the prior: although that selects `calc_ofv_ls()`, the optimizer data is
still built from `omega$est / weight_prior_var`, which is non-finite at zero
and fails in `solve()`. Use `type = "ls"` instead.

### Estimation Pipeline

```
get_map_estimates()
  ├── parse_weight_prior()      # -> weight_prior_var; picks calc_ofv_map vs calc_ofv_ls
  ├── check_inputs()            # validate all inputs upfront
  ├── parse_*()                 # ~10 parse_ functions normalize inputs
  │   ├── parse_input_data()
  │   ├── parse_omega_matrix()
  │   ├── parse_error()
  │   └── ...
  ├── mle_wrapper()             # wraps optim() + optional numDeriv Hessian
  │   └── ll_func_PKPDsim()     # likelihood: calls PKPDsim, applies error model + censoring
  ├── calc_residuals()          # g.o.f., CWRES, and the FOCE Jacobian
  ├── get_varcov_matrix()       # vcov, with omega as fallback
  └── get_mahalanobis()
```

Two details of this flow are non-obvious and easy to break:

- **Model class dispatch.** If the model does not carry a `cpp` attribute,
  `ll_func` silently falls back to `ll_func_generic()` with a warning instead
  of `ll_func_PKPDsim()`.
- **Hessian vs FOCE Jacobian.** When `residuals = TRUE`, `calc_residuals()`
  already computes a FOCE Jacobian that yields both CWRES and a vcov, so the
  more expensive `numDeriv::hessian()` in `mle_wrapper()` is skipped
  (`skip_hessian_mle <- skip_hessian || residuals`). The FOCE vcov
  (`obj$foce_vcov`) is then *preferred* over `obj$fit$vcov`. Changing either
  path affects both the reported vcov and CWRES.

### Mixture models

If the model carries a `mixture` attribute, `get_map_estimates()` fits the
model once per mixture value, converts the OFVs to posterior probabilities,
selects the most likely group, and returns the details under `obj$mixture`.
The population parameter is overwritten with the selected value.

### Key Function Roles

| File | Purpose |
|------|---------|
| `R/get_map_estimates.R` | Main user-facing function; orchestrates the full estimation pipeline |
| `R/mle_wrapper.R` | Wraps `optim()`, handles Hessian for variance-covariance |
| `R/ll_func_PKPDsim.R` | Computes log-likelihood by calling PKPDsim and applying residual error + censoring |
| `R/ll_func_generic.R` | Fallback likelihood for non-PKPDsim (non-`cpp`) models |
| `R/calc_ofv_map.R` | Objective function value for MAP (includes prior penalty) |
| `R/calc_ofv_ls.R` | OFV for least squares |
| `R/calc_residuals.R` | Residuals, g.o.f. metrics, and the FOCE Jacobian |
| `R/calc_cwres.R` | Conditional weighted residuals (FOCE approximation) |
| `R/get_varcov_matrix.R` | Builds the vcov output, falling back to omega |
| `R/parse_omega_matrix.R` | Converts various omega input formats to full covariance matrix |
| `R/parse_weight_prior.R` | Converts `weight_prior` (SD scale) to a variance scaling factor |
| `R/check_inputs.R` | Comprehensive upfront validation (throws descriptive errors) |
| `R/run_sequential_map.R` | Sequential MAP fits over time windows, for parameter "tracking" |
| `R/print.map_estimates.R` | Print method for the returned `map_estimates` object |

### Parameter Conventions

- **ETA (η)**: Random effects (inter-individual variability), parameterized on
  log-normal or normal scale. The `as_eta` argument marks parameters estimated
  directly on the eta scale.
- **Omega (Ω)**: Between-subject variability covariance matrix.
  `parse_omega_matrix()` accepts either a full matrix or a lower-triangle
  vector of variances/covariances, and also returns `$nonfixed` and starting
  `$eta`. CV% is not an accepted input: `create_block_from_cv()` is a separate
  helper that builds a diagonal lower-triangle block from a CV fraction.
- **IOV**: Inter-occasion variability, handled via `create_iov_object()` and
  the `iov_bins` argument.
- **Fixed parameters**: Listed in `fixed` argument to exclude from optimization.

### Return object

`get_map_estimates()` returns a list of class `map_estimates` containing
`fit`, `parameters` (individual estimates), `mixture`, `vcov_full` / `vcov`,
`mahalanobis`, and `prior` (the population parameters, omega, and fixed list
used). The residual/g.o.f. fields are added by `calc_residuals()` only when
`residuals = TRUE`; without them `mahalanobis` is `NULL`.

Only the non-mixture `mle_wrapper()` call is wrapped in a `tryCatch` that
returns the caught **error object itself** instead of throwing, so callers
should check the returned class. Mixture fits and errors raised elsewhere in
the pipeline propagate normally.

### Test Organization

Tests live in `tests/testthat/` and use testthat edition 3.

- `test-get_map_estimates.R` is the most comprehensive file; there are also
  focused variants (`.lagtime.R`, `.ss.R`) for lagtime and steady-state paths.
- Shared model fixtures (PKPDsim models built with `new_ode_model()`) are in
  `tests/testthat/setup.R` and are available to all test files automatically.
- `tests/testthat/output/` holds `verify_output()` regression files (e.g. the
  print-method output), referenced via `test_path("output", ...)`.
- `tests/testthat/nm/` contains NONMEM fixtures. `test-nonmem-vcov-comparison.R`
  cross-checks the FOCE vcov against NONMEM 7.5.1 using **hard-coded**
  reference values — NONMEM is licensed and is not run in CI. The control
  stream that produced those numbers is documented in the test header; update
  the header alongside any reference value.

## CI

GitHub Actions (`.github/workflows/R-CMD-check.yaml`) runs
`R CMD check --no-manual --as-cran` on Ubuntu/R-release on push/PR to `master`,
erroring only on `error` (not warnings). macOS, Windows, and R-devel runners are
present but commented out in the matrix. The workflow installs
`InsightRX/PKPDsim` as an extra package and needs `secrets.PAT_TOKEN` as
`GITHUB_PAT`.
