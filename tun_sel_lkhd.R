tun_sel_lkhd = function(pushed_A_est_by_lam, pushed_Beta_est_by_lam, 
                        data_info, ordin_pred_as = "num",
                        thre_list = c(0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 1)){
  
  if(F){ # Debugging input
    pushed_A_est_list = pushed_A_est_by_lam
    ordin_pred_as = "num"
    thre_list = c(0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 1)
  } 
  
  if(is.null(thre_list)) thre_list =  c(0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 1)
  
  types_by_node = data_info$types_by_node
  X1 = data_info$X1
  n_obs = data_info$n_obs
  X = cbind(X1, rep(1,n_obs))
  Y = data_info$Y
  rev_interv_info = data_info$rev_interv_info
  n_obs_by_node = apply(rev_interv_info,2,sum)
  rev_interv_info_design = data_info$rev_interv_info_design
  Sr = data_info$Sr
  T_H = data_info$T_H
  
  ## Only for mix.noteargis =====================================
  idx_resps_by_node = data_info$idx_resps_by_node 
  
  idx_conti_nodes = which(types_by_node == "c")
  n_obs_by_node_conti = n_obs_by_node[idx_conti_nodes]
  idx_multi_nodes = which(types_by_node == "m")
  idx_ordin_nodes = which(types_by_node == "o")
  
  idx_resps_multi = unlist(idx_resps_by_node[idx_multi_nodes])
  idx_resps_conti = unlist(idx_resps_by_node[idx_conti_nodes])
  idx_resps_ordin = unlist(idx_resps_by_node[idx_ordin_nodes])
  # =============================================================
  
  if (all(types_by_node == "o")) {
    node_type = "o"
  } else if (all(types_by_node == "m")) {
    node_type = "m"
  } else if (all(types_by_node == "c")) {
    node_type = "c"
  } else {
    node_type = "mix"
  }
  
  lkhd_list = 
    switch(node_type,
           "o" = {
             sapply(pushed_Beta_est_by_lam, function(x) {
               lkhd_ordin(x, X, Y, rev_interv_info, rev_interv_info_design, Sr, T_H)
             }) 
           },
           "m" = {
             sapply(pushed_Beta_est_by_lam, function(x) {
               lkhd_multi(x, X, Y, rev_interv_info, rev_interv_info_design, Sr)
             }) 
           },
           "c" = {
             sapply(pushed_Beta_est_by_lam, function(x) {
               lkhd_conti(x, X, Y, rev_interv_info)
             }) 
             
           },
           "mix" = {
             sapply(pushed_Beta_est_by_lam, function(x) {
               lkhd_mix(x, X, Y, rev_interv_info, rev_interv_info_design,
                        n_obs_by_node_conti, n_obs,
                        as.numeric(idx_resps_conti), as.numeric(idx_resps_multi), 
                        as.numeric(idx_resps_ordin),
                        as.numeric(idx_multi_nodes), as.numeric(idx_ordin_nodes),
                        Sr, T_H)
             })
           })
  
  n_edges = attr(pushed_A_est_by_lam, "n_edges")
  n_ests = length(pushed_A_est_by_lam)
  
  temp_list_set = cbind(n_edges, lkhd_list, 1:n_ests)
  sort_idx = order(n_edges, -lkhd_list)
  temp_list_sorted = temp_list_set[sort_idx,]
  temp_obj = temp_list_sorted[!duplicated(temp_list_sorted[,1]),]
  
  dr_list = abs(diff(temp_obj[,2]) /  diff(temp_obj[,1]))
  
  eval_K <- function(dr_list, thre){
    max(which(dr_list >= thre * max(dr_list)) + 1)
  }
  
  eval_K_vec <- Vectorize(eval_K,"thre")
  K_list_temp = eval_K_vec(dr_list, thre_list)
  K_list = temp_obj[K_list_temp, 3]
  names(K_list) = paste0("thre=",thre_list)
  
  return(K_list)
  
}

################################################################################
## Debuggig Block ##############################################################
################################################################################

if (F){
  Beta = pushed_Beta_est_by_lam[[3]]
  n_obs = nrow(X)
  pre_H = X%*%Beta
  H = pre_H %*% T_H
  exp_H = exp(H)
  log(1+ exp_H %*% t(Sr)) * rev_interv_info
  (Y*H)%*%t(Sr)
}

################################################################################
################################################################################
################################################################################