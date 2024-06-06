[![Build Status](https://magnum.travis-ci.com/InsightRX/PKPDmap.svg?token=qfpEFBKKaHdzzMBZjxnk&branch=master)](https://magnum.travis-ci.com/InsightRX/PKPDmap)

# PKPDmap

MAP Bayesian and non-parametric data fitting for PK(PD) models.

## Description

This R library implements various methods for fitting of individual PK(PD) data. For simulation of PK(PD) profiles it uses the `PKPDsim` library. The following estimation methods are implemented and available in the `get_map_estimates` function:

- `map`: MAP Bayesian: standard maximum a priori estimation of empirical Bayes estimates
- `map_flat_prior`: MAP Bayesian with flat prior: reduced weight of prior distribution, i.e. will put more trust in the data and hence reduce shrinkage (user can control the weighting).
- `ls`: Least Squares: not *truly* least squares but rather special case of MAP Bayesian with flat priors, i.e. with (almost) fully flat prior.
- `np`: non-parametric estimation. Performs non-parametric estimation based on user-specified parameter matrix of previously studied subjects. (Not exactly the same algorithm as in BestDose)
- `np_hybrid`: hybrid MAP and non-parametric extended grid. Will first perform standard MAP estimation (potentially with flattened prior), and subsequently perform non-parametric estimation on grid of support points around MAP estimates. (Not exactly the same algorithm as in BestDose)

## Sequential estimation

The module also supports sequential MAP Bayesian estimation (using anyone of the aforementioned methods), thus allowing the Bayes estimates to change over time. This can be useful for fitting data from patients in which PK is changing over time, such as commonly observed in patients on the ICU. This functionality is experimental and not fully developed yet.

## License

MIT
