# PKPDmap

MAP Bayesian fitting for PK(PD) models.

## Description

This R library implements various methods for fitting of individual PK(PD) data. For simulation of PK(PD) profiles it uses the `PKPDsim` [package](https://github.com/insightrx/pkpdsim). The following estimation methods are implemented and available in the `get_map_estimates` function via the `type` argument:

- `map`: MAP Bayesian: standard maximum a priori estimation of empirical Bayes estimates. For MAP estimation with flattened priors, please use the `weight_prior` argument.
- `pls`: An approximate implementation of penalized least squares, or more precisely, MAP Bayesian with a flattened prior. The prior weight is fixed internally to a very small value (regardless of the `weight_prior` argument).
- `ls`: Least squares: unlike `pls`, this is *not* a flat-prior special case of MAP. The objective function has no prior/omega term at all, and the residual error model is forced to `add = 1, prop = 0` regardless of the `error` argument, so it minimizes unweighted squared residuals.

## Installation

The development version of PKPDmap always has the most up-to-date improvements
and bug fixes. We aim to release PKPDmap on CRAN soon as well.

The development version of PKPDsim can be installed using:

```
devtools::install_github("InsightRX/PKPDmap")
```

## Contributing

We welcome input from the community:

- If you think you have encountered a bug, please [submit an issue](https://github.com/InsightRX/PKPDmap/issues) 
on the GitHub page. Please include a reproducible example of the unexpected 
behavior.

- Please [open a pull request](https://github.com/InsightRX/PKPDmap/pulls) if
you have a fix or updates that would improve the package. If you're not sure if
your proposed changes are useful or within scope of the package, feel free to
contact one of the authors of this package.

## Disclaimer

The functionality in this R package is provided "as is". While its authors 
adhere to software development best practices, the software may still contain 
unintended errors.

InsightRX Inc. and the authors of this package can not be held liable for any
damages resulting from any use of this software. By the use of this software 
package, the user waives all warranties, expressed or implied, including any 
warranties to the accuracy, quality or suitability of InsightRX for any 
particular purpose, either medical or non-medical.

<div align="right">
© <img src="man/figures/insightrx_logo_color.png" alt="InsightRX logo" width="120" />
</div>

