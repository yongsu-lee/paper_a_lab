#include "trexpm.h"
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
double lkhd_conti(mat Beta, mat X, mat Y, mat rev_interv_info) {
  
  // vec w0_temp = w.tail(n_resps); 
  // vec w1_temp = w.head(n_coefs);
  // mat W0_temp = reshape(w0_temp, 1, n_resps); 
  // mat W1_temp = reshape(w1_temp, idx_intcpt-1, n_resps);
  // mat W_temp = join_cols(W1_temp, W0_temp);
  // mat W = W_temp % white_coefs;
  // mat W1 = W.head_rows(idx_intcpt-1);
  
  double n_obs = X.n_rows;
  mat H = (X * Beta) % rev_interv_info;

  double max_lkhd = - (0.5) * accu(log(sum(pow(Y - H, 2),0)));
  
  return(max_lkhd);
  
}

