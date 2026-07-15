# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Package Overview

PKPDmap is an R package implementing Maximum A Posteriori (MAP) Bayesian estimation for pharmacokinetic/pharmacodynamic (PK/PD) data. It depends heavily on [PKPDsim](https://github.com/InsightRX/PKPDsim) (also an InsightRX package) for ODE-based PK/PD simulation.

## Common Commands

All development uses standard R/devtools workflows (no Makefile):

```r
# Load package for interactive development
devtools::load_all()

# Run all tests
devtools::test()

# Run a single test file
devtools::test(filter = "get_map_estimates")

# Regenerate documentation from roxygen2 comments
devtools::document()

# Full R CMD check (as done in CI)
devtools::check("--no-manual --as-cran")

# Install PKPDsim dependency from GitHub (requires PAT_TOKEN)
remotes::install_github("InsightRX/PKPDsim")
```

## Architecture

### Estimation Methods

The main entry point is `get_map_estimates()` (`R/get_map_estimates.R`), which supports five estimation methods via the `method` argument:

- `map` — standard MAP Bayesian with empirical Bayes (default)
- `map_flat_prior` — MAP with flattened priors (reduced shrinkage)
- `ls` — Least squares (nearly-flat MAP priors)


### Estimation Pipeline

```
get_map_estimates()
  ├── check_inputs()            # validate all inputs upfront
  ├── parse_*()/               # ~10 parse_ functions normalize inputs
  │   ├── parse_input_data()
  │   ├── parse_omega_matrix()
  │   ├── parse_error()
  │   └── ...
  └── mle_wrapper()            # wraps optim() + Hessian calculation
      └── ll_func_PKPDsim()   # likelihood: calls PKPDsim, applies error model + censoring
```

### Key Function Roles

| File | Purpose |
|------|---------|
| `R/get_map_estimates.R` | Main user-facing function; orchestrates the full estimation pipeline |
| `R/mle_wrapper.R` | Wraps `optim()`, handles Hessian for variance-covariance |
| `R/ll_func_PKPDsim.R` | Computes log-likelihood by calling PKPDsim and applying residual error + censoring |
| `R/calc_ofv_map.R` | Objective function value for MAP (includes prior penalty) |
| `R/calc_ofv_ls.R` | OFV for least squares |
| `R/parse_omega_matrix.R` | Converts various omega input formats to full covariance matrix |
| `R/check_inputs.R` | Comprehensive upfront validation (throws descriptive errors) |
| `R/run_sequential_map.R` | Batch MAP estimation across multiple individuals |

### Parameter Conventions

- **ETA (η)**: Random effects (inter-individual variability), parameterized on log-normal or normal scale
- **Omega (Ω)**: Between-subject variability covariance matrix — can be supplied as full matrix, CV%, or variance vector
- **IOV**: Inter-occasion variability, handled via `create_iov_object()`
- **Fixed parameters**: Listed in `fixed` argument to exclude from optimization

### Test Organization

Tests live in `tests/testthat/`. The most comprehensive test file is `test-get_map_estimates.R`. Snapshot outputs for regression testing are stored in `tests/testthat/output/`. Test setup fixtures (shared model objects) are in `tests/testthat/setup.R`.

## CI

GitHub Actions runs `R CMD check --no-manual --as-cran` on Ubuntu/R-release on push/PR to `master`. macOS and Windows runners are currently disabled in the matrix.
