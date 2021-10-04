#include "trexpm.h"
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
vec grad0cpp_mix(vec w, vec beta, vec mu, double rho, double alpha, double rho1,
                mat X, mat Y, mat white_coefs, 
                mat rev_interv_info,
                mat rev_interv_info_design, rowvec n_obs_by_node_conti,
                int idx_intcpt, int n_resps, 
                double n_obs, double n_nodes, int n_coefs, 
                vec idx_resps_conti, vec idx_resps_multi, vec idx_resps_ordin,
                vec idx_multi_nodes, vec idx_ordin_nodes,
                mat Sr, mat Se, mat SSr, vec rj_lev, vec rj_lev_sub,
                mat T_H, mat T_Mu, mat Y_tilde) {
  
  uvec cpp_rj_lev = conv_to<uvec>::from(rj_lev - 1);
  uvec cpp_rj_lev_sub = conv_to<uvec>::from(rj_lev_sub - 1);
  uvec cpp_ordin_resp_idx = conv_to<uvec>::from(idx_resps_ordin - 1);
  uvec cpp_multi_resp_idx = conv_to<uvec>::from(idx_resps_multi - 1);
  uvec cpp_conti_resp_idx = conv_to<uvec>::from(idx_resps_conti - 1);
  uvec cpp_multi_nodes = conv_to<uvec>::from(idx_multi_nodes - 1);
  uvec cpp_ordin_nodes = conv_to<uvec>::from(idx_ordin_nodes - 1);
  
  vec w0_temp = w.tail(n_resps); vec w1_temp = w.head(n_coefs);
  mat W0_temp = reshape(w0_temp, 1, n_resps); 
  mat W1_temp = reshape(w1_temp, idx_intcpt-1, n_resps);
  mat W_temp = join_cols(W1_temp, W0_temp);
  mat W = W_temp % white_coefs;
  mat W1 = W.head_rows(idx_intcpt-1);

  mat pre_H = X * W;
  mat H = pre_H % rev_interv_info_design;

  // grad for continuous part
  
  mat Y_conti = Y.cols(cpp_conti_resp_idx); 
  mat H_conti = H.cols(cpp_conti_resp_idx); 
  mat W_conti = W.cols(cpp_conti_resp_idx);

  rowvec obs_ratio = n_obs_by_node_conti / n_obs;
  rowvec Deno = sum(pow(Y_conti - H_conti, 2),0);
  mat grad_loss_conti_mat = - trans(X) * (Y_conti - H_conti);
  grad_loss_conti_mat.each_row() /= Deno;
  grad_loss_conti_mat.each_row() %= obs_ratio;

  // grad for multinomial part
  
  mat rev_interv_info_design_multi = rev_interv_info_design.cols(cpp_multi_resp_idx);
  mat Y_multi = Y.cols(cpp_multi_resp_idx);
  mat H_multi = H.cols(cpp_multi_resp_idx);
  mat W_multi = W.cols(cpp_multi_resp_idx);

  mat Mu_multi = exp(H_multi);
  mat Deno_P_multi = 1.0 + Mu_multi * trans(SSr);
  mat P_multi = (Mu_multi / Deno_P_multi) % rev_interv_info_design_multi;
  mat grad_loss_multi_mat =  -(1.0/n_obs) * trans(X) * (Y_multi - P_multi);
  
  // grad for ordinal part
  
  mat rev_interv_info_design_ordin = rev_interv_info_design.cols(cpp_ordin_resp_idx);
  mat Y_ordin = Y.cols(cpp_ordin_resp_idx);
  mat Y_tilde_ordin = Y_tilde.cols(cpp_ordin_resp_idx);
  mat pre_H_ordin = pre_H.cols(cpp_ordin_resp_idx);
  mat H_ordin = (pre_H_ordin * T_H) % rev_interv_info_design.cols(cpp_ordin_resp_idx);
  mat W_ordin = W.cols(cpp_ordin_resp_idx);
  mat exp_H_ordin = exp(H_ordin);
  
  mat Mu_tilde = exp_H_ordin * T_Mu;
  mat Deno_P_ordin = 1.0 + Mu_tilde.cols(cpp_rj_lev_sub);
  mat P_ordin = (Mu_tilde / Deno_P_ordin) % rev_interv_info_design_ordin;
  mat grad_loss_ordin_mat = -(1.0/n_obs) * trans(X) * (Y_tilde_ordin - P_ordin);
  
  // combined grad loss
  mat grad_loss_mix_mat(idx_intcpt, n_resps, fill::zeros);
  grad_loss_mix_mat.cols(cpp_conti_resp_idx) = grad_loss_conti_mat;
  grad_loss_mix_mat.cols(cpp_multi_resp_idx) = grad_loss_multi_mat;
  grad_loss_mix_mat.cols(cpp_ordin_resp_idx) = grad_loss_ordin_mat;
  mat grad_loss_mat = grad_loss_mix_mat % white_coefs; 

  // penalty part
  mat pp = Se * (W1 % W1) * trans(Sr);
  mat E = expmat(pp);
  double h = trace(E) - n_nodes;

  // grad_h_mat
  mat grad_h_mat(idx_intcpt, n_resps, fill::zeros);
  mat SetEtSr = trans(Se) * trans(E) * Sr;
  mat grad_h_mat_coef = SetEtSr % (2.0 * W1);
  grad_h_mat.head_rows(idx_intcpt-1) = grad_h_mat_coef;

  // grad_mat -> grad_vec
  mat grad_mat = grad_loss_mat + grad_h_mat * (rho1 * h + alpha);
  mat grad_mat1 = grad_mat.head_rows(idx_intcpt-1), grad_mat0 = grad_mat.tail_rows(1);
  vec grad_vec = join_cols(vectorise(grad_mat1), vectorise(grad_mat0));

  vec result = grad_vec + rho * (w - beta + mu);
  
  return(result);
  
}
