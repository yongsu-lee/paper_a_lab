#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
double arma_trexpm(mat A) {
  cx_vec eigval;
  eig_gen(eigval, A);
  cx_vec exp_eig = exp(eigval);
  double approx = real(sum(exp_eig));
  return(approx);
}

// [[Rcpp::export]]
mat arma_pp(mat Sr, mat Se, mat W1_new) {
  mat pp = Se * (W1_new % W1_new) * trans(Sr);
  return(pp);
}

// [[Rcpp::export]]
double arma_h_cal(vec w, mat white_coefs, int idx_intcpt, int n_resps, 
                  int n_coefs, mat Sr, mat Se, double n_nodes) {
  
  vec w0_temp = w.tail(n_resps); vec w1_temp = w.head(n_coefs);
  mat W0_temp = reshape(w0_temp, 1, n_resps); 
  mat W1_temp = reshape(w1_temp, idx_intcpt-1, n_resps);
  mat W_temp = join_cols(W1_temp, W0_temp);
  mat W = W_temp % white_coefs;
  mat W1 = W.head_rows(idx_intcpt-1);
  
  mat pp = Se * (W1 % W1) * trans(Sr);
  cx_vec eigval;
  eig_gen(eigval, pp);
  cx_vec exp_eig = exp(eigval);
  double approx = real(sum(exp_eig));
  double h_val = approx - n_nodes;
  return(h_val);
}

// W_new_temp = matricizing(w_new_temp, n_expls, n_resps)
//   W_new = W_new_temp * white_coefs
//   
//   W1_new = W_new[-idx_intcpt, ] # coef part only