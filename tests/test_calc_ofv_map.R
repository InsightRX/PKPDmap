library(testit)
library(PKPDmap)

assert(
  "Expected OFV",{
  eta <- rep(0, 5)
  omega <- structure(
    c(
      0.284, 0, 0, 0, 0, 
      0, 0.0868, 0.0464, 0.0799, 0.11, 
      0, 0.0464, 0.0339, 0.0511, 0.0953, 
      0, 0.0799, 0.0511, 0.524, 0.395, 
      0, 0.11, 0.0953, 0.395, 0.44
    ), 
    .Dim = c(5L, 5L)
  )
  dv <- c(4.6, 3.9)
  ipred <- c(3.3, 3.1)
  res_sd <- rep(0.29, 2)
  
  out <- calc_ofv_map(eta, omega, dv, ipred, res_sd)
  all(round(out, 4) == c(3.3149, -9.7286, -3.4861))
})
