gen_lambdas = function(data_info, n_lams = 30, eps = 0.05, ver = "1"){
  
  if (F){
    n_lams = 5; eps = 0.05; ver = "1"
  }
  
  Sr = data_info$Sr
  Se = data_info$Se
  n_obs = data_info$n_obs
  X1 = data_info$X1
  Y = data_info$Y
  rev_interv_info_design = data_info$rev_interv_info_design
  n_obs_by_node = apply(rev_interv_info_design,2,sum)

  types_by_node = data_info$types_by_node
  conti_nodes = which(types_by_node == "c")
  multi_nodes = which(types_by_node == "m")
  ordin_nodes = which(types_by_node == "o")
  idx_resps_by_node = data_info$idx_resps_by_node
  #n_idx_expls_by_node = data_info$idx_expls_by_node

  Ybar = apply(Y, 2, sum) / n_obs_by_node
  Yc_temp = sweep(Y, 2, Ybar, "-")
  Yc = Yc_temp * rev_interv_info_design
  X1tYc = t(X1) %*% Yc

  ## continuous part ####
  lambda_max_conti = NULL
  
  if (any(conti_nodes)) {
    
    conti_resp = unlist(idx_resps_by_node[(types_by_node == "c")])
    
    Yc_norm_sq = apply(Yc, 2, function(x) sum(x^2))
    scaled_X1tYc = sweep(X1tYc, 2, Yc_norm_sq, "/" )
    norm_scaled_X1tYc = sqrt(Se %*% (scaled_X1tYc^2) %*% t(Sr))
    diag(norm_scaled_X1tYc) <- 0

    lambda_max_conti = max(norm_scaled_X1tYc[,conti_nodes, drop = F])
  }
  
  ## multinomial part ####
  ### Note: ordinal part can be merged here.
  lambda_max_multi = NULL
  
  if (any(multi_nodes)){
    
    norm_X1tYc = sqrt(Se %*% (X1tYc^2) %*% t(Sr))
    diag(norm_X1tYc) <- 0
    
    multi_resp = unlist(idx_resps_by_node[multi_nodes])
    lambda_max_multi = max( (1/n_obs) * norm_X1tYc[, multi_nodes, drop = F] )
    
  }
  
  ## ordinal part ####
  lambda_max_ordin = NULL
  
  if (any(ordin_nodes)){

    Y_tilde = data_info$Y_tilde
    Y_tilde_bar = apply(Y_tilde, 2, sum) / n_obs_by_node
    Y_tilde_c_temp = sweep(Y_tilde, 2, Y_tilde_bar, "-")
    Y_tilde_c = Y_tilde_c_temp * rev_interv_info_design
    X1tY_tilde_c = t(X1) %*% Y_tilde_c 
    
    norm_X1tY_tilde_c = sqrt(Se %*% (X1tY_tilde_c^2) %*% t(Sr))
    diag(norm_X1tY_tilde_c) <- 0
    
    # ordin_resp = unlist(idx_resps_by_node[ordin_nodes]) ## ?????
    lambda_max_ordin = max( (1/n_obs) * norm_X1tY_tilde_c[,ordin_nodes, drop = F] )
    
  }
  
  
  ## find the max ####
  lambda_max = max(lambda_max_conti, lambda_max_multi, lambda_max_ordin) * 1.01 
  # ... *1.01 adjust lambda_max slightly larger to make it enough larger all the
  # ... all the coefs are zeros
  
  ## generate list of lambdas ####
  lambda_min = lambda_max * eps
  
  if (ver == "1"){
    lambdas = exp(seq(log(lambda_max), log(lambda_min), length.out = n_lams))
  } else if (ver == "2") {
    lambdas = seq(lambda_max, lambda_min, length.out = n_lams)
  }
  
  return(lambdas)
  
}
