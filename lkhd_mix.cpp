#include "trexpm.h"
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
double lkhd_mix(mat Beta, mat X, mat Y, 
                mat rev_interv_info, mat rev_interv_info_design, 
                rowvec n_obs_by_node_conti, double n_obs,
                vec idx_resps_conti, vec idx_resps_multi, vec idx_resps_ordin,
                vec idx_multi_nodes, vec idx_ordin_nodes,
                mat Sr, mat T_H) {
  
  mat pre_H = X * Beta;
  mat H = pre_H % rev_interv_info_design;
  
  uvec cpp_conti_resp_idx = conv_to<uvec>::from(idx_resps_conti - 1);
  uvec cpp_multi_resp_idx = conv_to<uvec>::from(idx_resps_multi - 1);
  uvec cpp_ordin_resp_idx = conv_to<uvec>::from(idx_resps_ordin - 1);
  
  uvec cpp_multi_nodes = conv_to<uvec>::from(idx_multi_nodes - 1);
  uvec cpp_ordin_nodes = conv_to<uvec>::from(idx_ordin_nodes - 1);
  
  mat Y_conti = Y.cols(cpp_conti_resp_idx); 
  mat Y_multi = Y.cols(cpp_multi_resp_idx); 
  mat Y_ordin = Y.cols(cpp_ordin_resp_idx); 
  mat H_conti = H.cols(cpp_conti_resp_idx); 
  mat H_multi = H.cols(cpp_multi_resp_idx); 
  
  // Continuous case
  
  rowvec loglkhds_conti = log(sum(pow(Y_conti - H_conti, 2),0));
  double loss_conti = - (0.5) * (1.0/n_obs) * accu(loglkhds_conti % n_obs_by_node_conti);
  
  // Multinomial case
  
  mat Sr_multi = Sr.submat(cpp_multi_nodes, cpp_multi_resp_idx);
  mat rev_interv_info_multi = rev_interv_info.cols(cpp_multi_nodes);
  mat Mu_multi = exp(H_multi);
  mat Deno_multi = log(1.0 + Mu_multi * trans(Sr_multi)) % rev_interv_info_multi;
  mat Nume_multi = (Y_multi % H_multi) * trans(Sr_multi);
  double loss_multi = (1.0/n_obs) * accu(Nume_multi - Deno_multi);
  
  // Ordinal case
  
  mat pre_H_ordin = pre_H.cols(cpp_ordin_resp_idx);
  mat H_ordin = (pre_H_ordin * T_H) % rev_interv_info_design.cols(cpp_ordin_resp_idx);
  mat exp_H_ordin = exp(H_ordin);
  
  mat Sr_ordin = Sr.submat(cpp_ordin_nodes, cpp_ordin_resp_idx);
  mat rev_interv_info_ordin = rev_interv_info.cols(cpp_ordin_nodes);
  mat Deno_ordin = log(1.0 + exp_H_ordin * trans(Sr_ordin)) % rev_interv_info_ordin;
  mat Nume_ordin = (Y_ordin % H_ordin) * trans(Sr_ordin);
  double loss_ordin = (1.0/n_obs) * accu(Nume_ordin - Deno_ordin);

  double max_lkhd = (loss_conti + loss_multi + loss_ordin); 
  
  return(max_lkhd);
  
}