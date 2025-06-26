# Changes in master compared to last release (20240610)

- many code chunks in main estimation function now separated off to their own function
- improved handling of floating points during warning check
- The print method for map_estimates() now uses 
  `sqrt(diag(PKPDsim::triangle_to_full(x$prior$omega$full)))` instead of 
  `sqrt(diag(PKPDsim::triangle_to_full(x$prior$omega)))` for etas
- fix for mixture models: ipreds for mixture models were found not to be correct
- censoring information is now returned with the fit object
