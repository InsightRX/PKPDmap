// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec Mahalanobis(arma::mat x, arma::rowvec center, arma::mat cov){
  int n = x.n_rows;
  arma::mat x_cen;
  x_cen.copy_size(x);
  for (int i=0; i < n; i++) {
    x_cen.row(i) = x.row(i) - center;
  }
  return sum((x_cen * cov.i()) % x_cen, 1);    
}

// [[Rcpp::export]]
Rcpp::NumericVector dmvnorm_arma(arma::vec x,  arma::rowvec mean,  arma::mat sigma, bool log = false) { 
  arma::mat y(x);
  y.reshape(1, x.size());
  arma::vec distval = Mahalanobis(y,  mean, sigma);
  double logdet = sum(arma::log(arma::eig_sym(sigma)));
  double log2pi = std::log(2.0 * M_PI);
  arma::vec logretval = -( (y.n_cols * log2pi + logdet + distval)/2  ) ;
  
  if (log){
    return(Rcpp::NumericVector(logretval.begin(), logretval.end()));
  } else {
    return(exp(Rcpp::NumericVector(logretval.begin(), logretval.end())));
  }
}
