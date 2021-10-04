glmdag = function(data_input, types_by_node = NULL, ordin_pred_as = "num", 
                  interv_info=NULL, root_nodes=NULL, terminal_nodes=NULL,
                  lambdas=NULL, n_lams=NULL, 
                  path_par = FALSE, path_par_num = NULL,
                  thre_list = NULL,
                  admm_args = admm_arg_ctrl(), opt_method = "lbfgs",
                  fac_grp_lasso=FALSE, add_stop_rule=TRUE,
                  fit_hist=FALSE, verbose=FALSE, 
                  debug_mode=FALSE){
  
  if(F){ # debugging block ####
    ordin_pred_as = "num";
    interv_info=NULL; root_nodes=NULL; terminal_nodes=NULL;
    lambdas=NULL; n_lams=3; thre_list = NULL;
    path_par = FALSE; path_par_num = NULL;
    admm_args = admm_arg_ctrl(); opt_method = "lbfgs";
    fac_grp_lasso=FALSE; add_stop_rule=TRUE; fit_hist=FALSE;
    verbose=T; debug_mode=F
  }
  
  if (!is.data.frame(data_input)){
    stop("Input data should be a dataframe.")
  }
  
  if (is.null(types_by_node)){
    message("'types_by_node' is not provided. 'types_by_node' will be determined based on str(data_input) information.")
    types_by_node = ifelse((sapply(data_input, class) %in% "factor"), "m","c")
    ord_node = which(sapply(data_input, function(x) any(class(x) %in% "ordered")))
    types_by_node[ord_node] <- "o"
  }
  
  ## Normalizing for continuous variables
  idx_conti_nodes = which(types_by_node == "c")
  temp_scale = scale(data_input[idx_conti_nodes])
  data_input[idx_conti_nodes] <- temp_scale
  attr(data_input, "conti_mean") <- attr(temp_scale, "scaled:center")
  attr(data_input, "conti_sd") <- attr(temp_scale, "scaled:scale")
  
  ## Check if the input is proper ####
  # data_input = check_input(data_input)
  
  ## Read and set data ####
  n_obs = nrow(data_input)
  n_nodes = ncol(data_input)
  
  if (is.null(interv_info)){
    interv_info = matrix(FALSE, n_obs, n_nodes)
  }
  
  n_levels_by_node_temp = sapply(data_input, nlevels)
  n_levels_by_node = ifelse(types_by_node == "c", 1, n_levels_by_node_temp)
  
  n_resps_by_node = gen_n_xxxxx_by_node("resps", types_by_node, n_levels_by_node)
  n_expls_by_node = gen_n_xxxxx_by_node("expls", types_by_node, n_levels_by_node)
  
  idx_resps_by_node = gen_indices_by_node(n_resps_by_node)
  idx_expls_by_node = gen_indices_by_node(n_expls_by_node)
  
  n_expls <- sum(n_expls_by_node)
  n_resps <- sum(n_resps_by_node)
  n_coefs <- n_expls * n_resps
  
  ## Generate design matrices
  X1      = gen_X1_design(data_input, types_by_node, n_expls_by_node)
  Y_temp  = gen_Y_design(data_input, types_by_node, n_resps_by_node )
  
  ## Parameter modifier
  Sr = t(diag(n_nodes)[rep(1:n_nodes, times = n_resps_by_node),])
  Se = t(diag(n_nodes)[rep(1:n_nodes, times = n_expls_by_node),])
  
  ## Parameter modifier - nominal data only
  n_multi_nodes = sum(types_by_node == "m")
  n_resps_multi_only =  n_resps_by_node[types_by_node == "m"]
  if (n_multi_nodes >= 1){
    SSr_temp = 
      t(diag(n_multi_nodes)[rep(1:n_multi_nodes, times =  n_resps_multi_only),])
    SSr = SSr_temp[rep(1:n_multi_nodes, times =  n_resps_multi_only),]
  } else {
    SSr = matrix(0,0,0)
  }
  
  ## Parameter modifier - Ordinal data only
  idx_resps_ordin_only = idx_resps_by_node[types_by_node == "o"]
  n_resps_ordin_only = n_resps_by_node[types_by_node == "o"]
  rj_lev = unlist(lapply(idx_resps_ordin_only, function(x) rep(x[length(x)],length(x))))
  temp_diag = as.list(n_resps_ordin_only)
  if (any(types_by_node == "o")){
    T_H = bdiag(lapply(temp_diag, function(x) 1*lower.tri(matrix(1,x,x), diag=T)))
    T_Mu = bdiag(lapply(temp_diag, function(x) 1*upper.tri(matrix(1,x,x), diag=T)))
    T_H <- as.matrix(T_H); T_Mu <- as.matrix(T_Mu)
  } else {
    T_H <- T_Mu <- matrix(nrow = 0, ncol = 0 )
  }
  
  # Note rj_lev_sub should be defined in the mixed case
  if (any(types_by_node == "o")){
    types_temp = rep("o", length(n_resps_ordin_only))
    n_levels_temp = n_levels_by_node[types_by_node == "o"]
    n_resps_temp = gen_n_xxxxx_by_node("resps", types_temp, n_levels_temp)
    idx_resps_temp = gen_indices_by_node(n_resps_temp)
    rj_lev_sub = unlist(lapply(idx_resps_temp, function(x) rep(x[length(x)],length(x))))
  } else {
    rj_lev_sub = numeric()
  }
  
  n_expls = sum(n_expls_by_node)
  n_resps = sum(n_resps_by_node)
  n_coefs = n_expls * n_resps
  
  ## Intervention data adjustment
  rev_interv_info= !interv_info
  rev_interv_info_design = rev_interv_info[,rep(1:n_nodes, times = n_resps_by_node)]
  Y = Y_temp * rev_interv_info_design
  
  ## Response modifier - ordin only
  Y_tilde = Y
  if (any(types_by_node == "o")){
    vec_idx_resps_ordin_only = unlist(idx_resps_ordin_only)
    Y_ordin_only = Y[, vec_idx_resps_ordin_only] %*% T_Mu 
    Y_tilde[, vec_idx_resps_ordin_only] <- Y_ordin_only
  }
  
  data_info = list( n_resps_by_node = n_resps_by_node, 
                    types_by_node = types_by_node,
                    X1 = X1, Y = Y,
                    idx_resps_by_node = idx_resps_by_node,
                    idx_expls_by_node = idx_expls_by_node ,
                    n_obs = n_obs, n_nodes = n_nodes, n_coefs = n_coefs,
                    n_expls = n_expls, n_resps = n_resps,
                    Sr = Sr, Se = Se, SSr = SSr, 
                    rj_lev = rj_lev, rj_lev_sub = rj_lev_sub,
                    rev_interv_info = rev_interv_info,
                    rev_interv_info_design = rev_interv_info_design,
                    T_H = T_H, T_Mu = T_Mu, Y_tilde = Y_tilde)
  
  ## Set tuning parameters ####
  if (is.null(lambdas)){
    if (is.null(n_lams)){
      n_lams = 30
      lambdas = gen_lambdas(data_info, n_lams = n_lams)
    } else {
      lambdas = gen_lambdas(data_info, n_lams = n_lams)
    }
  } else {
    if (is.null(n_lams)) {
      n_lams = length(lambdas)
    } else {
      if (length(lambdas) != n_lams)
        n_lams = length(lambdas)
    }
  }
  
  if (path_par == T){
    if (is.null(path_par_num)) stop("path_par = TRUE requires path_par_num.")
    lambdas = lambdas[path_par_num]
  }
  
  ## Set W_update according to types_by_node ####
  W_update =  mix.W_update
  
  ## This part is designed for the final step (slightly faster)
  # W_update = switch(data_type, "multi-level" = multi.W_update, 
  #                   "gaussian" = gauss.W_update,
  #                   "mix" = mix.W_update)
  
  ## Set penalizing method 
  Beta_update = grp_lasso.Beta_update
  ### ... Later, can be modified by allowing various penalizing
  
  
  ## Initialize parameters ####
  terminal_idx = which(names(data_input) %in% terminal_nodes)
  root_idx = which(names(data_input) %in% root_nodes)
  
  init_W_info = initialize_W(data_info, terminal_idx, root_idx)
  W_init = init_W_info$W_init; white_coefs = init_W_info$white_coefs;
  Beta_init = W_init;
  Mu = matrix(0.0, n_expls + 1, n_resps)
  
  
  ## Debugging mode output ####
  if (debug_mode == T) {
    W_update <<- W_update; Beta_update <<- Beta_update
    W_init <<- W_init; Beta_init <<- Beta_init; white_coefs <<- white_coefs
    data_info <<- data_info; lambdas <<- lambdas; admm_args <<- admm_args;
    add_stop_rule <<- add_stop_rule; fac_grp_lasso <<- fac_grp_lasso
    verbose <<- verbose; fit_hist <<- fit_hist;
    message("Debugging mode is activated. No actual fitting is conducted.")
    return(invisible(NULL))
  }
  
  
  ## ADMM loop ####
  admm_results = admm_loop(W_update, Beta_update, W_init, Beta_init, white_coefs,
                           data_info, lambdas, admm_args, 
                           add_stop_rule, fac_grp_lasso, 
                           verbose, fit_hist, debug_mode)
  
  ## Pushed DAG 
  modified_dag_result = push_dag(admm_results$Beta_est_by_lam, data_info)
  
  if (path_par == T){
    
    pushed_A_est_by_lam = modified_dag_result$A_est_by_lam
    pushed_Beta_est_by_lam = modified_dag_result$Beta_est_by_lam
    
    fit = c(list(graph_est_by_lam = modified_dag_result$graph_est_by_lam, 
                 pushed_A_est_by_lam = pushed_A_est_by_lam,
                 pushed_Beta_est_by_lam = pushed_Beta_est_by_lam), 
            admm_results, data_info)
    
    
  } else {
    
    ## Select tuning parameters
    pushed_A_est_by_lam = modified_dag_result$A_est_by_lam
    pushed_Beta_est_by_lam = modified_dag_result$Beta_est_by_lam
    sel_tun = tun_sel_lkhd(pushed_A_est_by_lam, pushed_Beta_est_by_lam, data_info)
    
    fit = c(list(graph_est_by_lam = modified_dag_result$graph_est_by_lam, 
                 sel_tun = sel_tun,
                 pushed_A_est_by_lam = pushed_A_est_by_lam,
                 pushed_Beta_est_by_lam = pushed_Beta_est_by_lam), 
            admm_results, data_info)
    
  }
  
  class(fit) <- "glmdag"
  return(fit)
  
}