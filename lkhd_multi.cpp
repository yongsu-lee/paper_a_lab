#include "trexpm.h"
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
double lkhd_multi(mat Beta, mat X, mat Y, 
                  mat rev_interv_info, mat rev_interv_info_design,
                  mat Sr) {
  
  // vec w0_temp = w.tail(n_resps); 
  // vec w1_temp = w.head(n_coefs);
  // mat W0_temp = reshape(w0_temp, 1, n_resps); 
  // mat W1_temp = reshape(w1_temp, idx_intcpt-1, n_resps);
  // mat W_temp = join_cols(W1_temp, W0_temp);
  // mat W = W_temp % white_coefs;
  // mat W1 = W.head_rows(idx_intcpt-1);
  
  double n_obs = X.n_rows;
  mat H = (X * Beta) % rev_interv_info_design;
  mat Mu_multi = exp(H);
  
  mat Deno = log(1.0 + Mu_multi * trans(Sr)) % rev_interv_info;
  mat Nume = (Y % H) * trans(Sr);
  
  double max_lkhd = (1.0/n_obs) * accu(Nume - Deno) ;
  
  return(max_lkhd);
  
}

