initialize_W = function(data_info, terminal_idx, root_idx){
  
  ## Load necessary data
  Y = data_info$Y
  types_by_node = data_info$types_by_node
  rev_interv_info_design = data_info$rev_interv_info_design
  n_obs_by_node = apply(rev_interv_info_design,2,sum)
  
  n_nodes = data_info$n_nodes
  idx_resps_by_node = data_info$idx_resps_by_node
  idx_expls_by_node = data_info$idx_expls_by_node 
  
  idx_intcpt = data_info$n_expls + 1 ## including intercepts 
  n_resps = data_info$n_resps
  
  ## Initialize containers
  W_init = matrix(0, ncol = n_resps, nrow = idx_intcpt)
  white_coefs = matrix(1, ncol = n_resps, nrow = idx_intcpt)
  
  ## Caching Y_mean
  temp_Y_mean = apply(Y, 2, sum) / n_obs_by_node
  
  ## Continuous nodes ####
  idx_resps_conti_only = unlist(idx_resps_by_node[(types_by_node == "c")]) 
  init_est_conti = temp_Y_mean[idx_resps_conti_only]
  W_init[idx_intcpt, idx_resps_conti_only] = init_est_conti  
  
  ## Multinomial nodes ####
  idx_resps_multi_only = unlist(idx_resps_by_node[(types_by_node == "m")])
  idx_multi_nodes = idx_resps_by_node[(types_by_node == "m")]
  temp_split = lapply(idx_multi_nodes, function(x) apply(Y, 2, mean)[x] )
  
  init_est_multi = unlist(lapply(temp_split, function(x) log(x/(1-sum(x)))))
  W_init[idx_intcpt, idx_resps_multi_only] = init_est_multi
  
  ## Ordinal nodes ####
  idx_resps_ordin_only = unlist(idx_resps_by_node[(types_by_node == "o")])
  idx_ordin_nodes = idx_resps_by_node[(types_by_node == "o")]
  temp_split = lapply(idx_ordin_nodes, function(x) apply(Y, 2, mean)[x] )
  
  init_est_ordin = unlist(lapply(temp_split, function(x) {
      temp = c(x, 1-sum(x))
      log(x / temp[-1]) } ))
  W_init[idx_intcpt, idx_resps_ordin_only] = init_est_ordin
  
  ## Set white coefficients
  for (j in 1:n_nodes){ # j = 1
    
    idx_expls_j = idx_expls_by_node[[j]];
    idx_resps_j = idx_expls_by_node[[j]];
    
    if (any(terminal_idx %in% j)){
      white_coefs[idx_expls_j, ] <- 0
      
    } else if (any(root_idx %in% j)){
      white_coefs[,idx_resps_j ] <- 0
      
    } else {
      white_coefs[idx_expls_j, idx_resps_j] <- 0
    }
    
  }
  
  init_W_info = list(W_init = W_init, white_coefs = white_coefs)
  
  return(init_W_info)
  
}