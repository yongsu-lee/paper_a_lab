#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
double trexpm(mat A) {
  cx_vec eigval;
  eig_gen(eigval, A);
  cx_vec exp_eig = exp(eigval);
  double approx = real(sum(exp_eig));
  
  return(approx);
}

