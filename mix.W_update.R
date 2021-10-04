mix.W_update = function(W, Beta, Mu, rho, data_info, white_coefs, verbose,
                        rho1_max = 1e+20, h_tol = 1e-8, 
                        debug_mode = FALSE, ver = "1", opt_method = "lbfgs"){
  
  if (F){
    rho1_max = 1e+20; h_tol = 1e-8; debug_mode = FALSE; 
    ver = "1"; opt_method = "lbfgs"
  }
  
  ## Read data ####
  n_expls = data_info$n_expls
  idx_intcpt = data_info$n_expls + 1
  n_resps = data_info$n_resps
  X1 = data_info$X1
  
  Y = data_info$Y
  Y_tilde = data_info$Y_tilde  # ordinal data only
  rev_interv_info = data_info$rev_interv_info
  rev_interv_info_design = (data_info$rev_interv_info_design) * 1
  n_obs_by_node = apply(rev_interv_info,2,sum)
 
  Sr = data_info$Sr
  Se = data_info$Se
  SSr = data_info$SSr
  rj_lev = data_info$rj_lev
  rj_lev_sub = data_info$rj_lev_sub
  T_H = data_info$T_H
  T_Mu = data_info$T_Mu
  
  n_obs = data_info$n_obs
  n_coefs = data_info$n_coefs
  n_nodes = data_info$n_nodes
  
  X = cbind(X1, rep(1,n_obs))
  colnames(X)[idx_intcpt] <- "(Intcpt)"
  
  ## Only for mix.noteargis =====================================
  n_resps_by_node = data_info$n_resps_by_node
  types_by_node = data_info$types_by_node
  
  idx_resps_by_node = data_info$idx_resps_by_node 
  
  idx_conti_nodes = which(types_by_node == "c")
  n_obs_by_node_conti = n_obs_by_node[idx_conti_nodes]
  idx_multi_nodes = which(types_by_node == "m")
  idx_ordin_nodes = which(types_by_node == "o")
  
  idx_resps_multi = unlist(idx_resps_by_node[idx_multi_nodes])
  idx_resps_conti = unlist(idx_resps_by_node[idx_conti_nodes])
  idx_resps_ordin = unlist(idx_resps_by_node[idx_ordin_nodes])
  # =============================================================
  
  ## Initialize parameters ####
  W = W * white_coefs;
  w = vectorizing(W, idx_intcpt)
  beta = vectorizing(Beta, idx_intcpt)
  mu = vectorizing(Mu, idx_intcpt)
  
  alpha = 0.0
  rho1 = 1.0
  h = Inf
  
  if (debug_mode ==T) { 
    w <<- w; assign("beta", beta, envir = .GlobalEnv); mu <<- mu; 
    alpha <<- alpha; rho1 <<- rho1
    X <<- X; Y <<- Y
    idx_intcpt <<- idx_intcpt; n_resps <<- n_resps; n_coefs <<- n_coefs
    n_obs <<- n_obs; Sr <<- Sr; Se <<- Se; SSr <<- SSr;
    rj_lev <<- rj_lev;rj_lev_sub <<- rj_lev_sub; T_H <<-T_H;  T_Mu <<-  T_Mu
    Y_tilde <<- Y_tilde
    
    idx_multi_nodes <<- idx_multi_nodes; idx_ordin_nodes <<- idx_ordin_nodes;
    n_obs_by_node_conti <<- n_obs_by_node_conti;
    
    idx_resps_conti <<- idx_resps_conti; idx_resps_multi <<- idx_resps_multi;
    idx_resps_ordin <<- idx_resps_ordin; 
    
    rev_interv_info_design <<- rev_interv_info_design
    rev_interv_info <<- rev_interv_info
    message("Debugging mode is activated. No actual fitting is conducted.")
    return(invisible(NULL))
  }
  
  ## Debugging obj cpp function ####
  
  if (F){
    set.seed(10)
    W_temp = matrix(rnorm(length(w), mean = 0.01, sd = 0.01), nrow(W), ncol(W))
    W = W_temp * white_coefs;
    
    w = vectorizing(W * white_coefs, idx_intcpt)
    W1 = W[-idx_intcpt,]
    pre_H = (X%*%W)
    H = pre_H * rev_interv_info_design
    
    Y_conti = Y[, idx_resps_conti, drop = F]
    H_conti = H[, idx_resps_conti, drop = F]
    
    loglikhds_conti = log(apply(Y_conti - H_conti, 2, function(x) sum(x^2)))
    loss_conti = 0.5 * (1/n_obs) * sum(loglikhds_conti * n_obs_by_node_conti)
    # (1/n_obs) is reserved for the interventional case later on
    
    Y_multi = Y[, idx_resps_multi, drop = F]
    H_multi = H[, idx_resps_multi, drop = F]
    
    Sr_multi = Sr[idx_multi_nodes, idx_resps_multi]
    Mu_multi = exp(H_multi)
    rev_interv_info_multi = rev_interv_info[, idx_multi_nodes, drop = F]
    Deno_multi = log(1 + Mu_multi %*% t(Sr_multi)) * rev_interv_info_multi
    Nume_multi = (Y_multi * H_multi) %*% t(Sr_multi)
    loss_multi = - (1/n_obs) * sum(Nume_multi-Deno_multi)
    
    Y_ordin = Y[, idx_resps_ordin, drop = F]
    pre_H_ordin = pre_H[, idx_resps_ordin, drop =F]
    H_ordin = (pre_H_ordin %*% T_H) * rev_interv_info_design[, idx_resps_ordin, drop = F]
    exp_H_ordin = exp(H_ordin)
    
    Sr_ordin = Sr[idx_ordin_nodes, idx_resps_ordin]
    rev_interv_info_ordin = rev_interv_info[, idx_ordin_nodes]
    
    Deno_ordin = log( 1 + exp_H_ordin %*% t(Sr_ordin)) * rev_interv_info_ordin
    Nume_ordin = (Y_ordin * H_ordin) %*% t(Sr_ordin)
    loss_ordin = - (1/n_obs) * sum(Nume_ordin - Deno_ordin)
    
    pp = Se %*% (W1^2) %*% t(Sr)
    h = sum(diag(expm(pp))) - n_nodes
    vector = w - beta + mu
    loss = loss_conti + loss_multi + loss_ordin
    loss + (rho1/2) * h^2 + alpha * h + (rho/2) * sum(vector^2)
    
    sourceCpp("./obj0_mix.cpp")
    
  }
  
  ## Wrapping obj function ####
  obj1 = function(w, alpha, rho1){
    obj0cpp_mix(w, beta, mu, rho, alpha, rho1, X, Y, white_coefs,
                rev_interv_info, rev_interv_info_design, n_obs_by_node_conti,
                idx_intcpt, n_resps, n_obs, n_nodes, n_coefs, 
                as.numeric(idx_resps_conti), as.numeric(idx_resps_multi), 
                as.numeric(idx_resps_ordin),
                as.numeric(idx_multi_nodes), as.numeric(idx_ordin_nodes),
                Sr, Se, T_H)}
  
  ## Debugging obj cpp function ####
  if (F){
    
    obj1(w, alpha, rho1)
    
    set.seed(10)
    W_temp = matrix(rnorm(length(w), mean = 0.01, sd = 0.01), nrow(W), ncol(W))
    W = W_temp * white_coefs;
    
    
    w = vectorizing(W * white_coefs, idx_intcpt)
    W1 = W[-idx_intcpt,]
    pre_H = (X%*%W)
    H = pre_H * rev_interv_info_design
    
    ## Continuous nodes
    Y_conti = Y[, idx_resps_conti, drop = F]
    H_conti = H[, idx_resps_conti, drop = F]
    W_conti = W[, idx_resps_conti, drop = F] 
    
    obs_ratio = n_obs_by_node_conti / n_obs
    Deno_conti = apply(Y_conti - H_conti, 2, function(x) sum(x^2))
    grad_loss_conti_mat_temp0 = - t(X) %*% (Y_conti - H_conti)
    grad_loss_conti_mat_temp1 = sweep(grad_loss_conti_mat_temp0, 2, Deno_conti, "/")
    grad_loss_conti_mat = sweep(grad_loss_conti_mat_temp1, 2, obs_ratio, "*")
    
    ## Multinomial nodes
    Y_multi = Y[, idx_resps_multi, drop = F]
    H_multi = H[, idx_resps_multi, drop = F]
    W_multi = W[, idx_resps_multi, drop = F]
    rev_interv_info_design_multi = rev_interv_info_design[, idx_resps_multi]
    
    Mu_multi = exp(H_multi)
    Deno_P_multi = 1 + Mu_multi %*% t(SSr)
    P_multi = (Mu_multi / Deno_P_multi) * rev_interv_info_design_multi
    grad_loss_multi_mat = -(1/n_obs) * t(X) %*% (Y_multi - P_multi)
    
    ## Ordinal nodes
    Y_ordin = Y[, idx_resps_ordin]
    Y_tilde_ordin = Y_tilde[, idx_resps_ordin, drop=F]
    pre_H_ordin = pre_H[, idx_resps_ordin, drop = F]
    H_ordin = (pre_H_ordin %*% T_H) * rev_interv_info_design[, idx_resps_ordin, drop = F]
    W_ordin = W[, idx_resps_ordin]
    rev_interv_info_design_ordin = rev_interv_info_design[, idx_resps_ordin, drop = F]
    exp_H_ordin = exp(H_ordin)
    
    Mu_tilde = exp_H_ordin %*% T_Mu
    Deno_P_ordin = 1 + Mu_tilde[, rj_lev_sub]
    P_ordin = (Mu_tilde / Deno_P_ordin) * rev_interv_info_design_ordin
    grad_loss_ordin_mat = -(1/n_obs) * t(X) %*% (Y_tilde_ordin - P_ordin)

    ## Combine grad_loss_#####_mat 
    grad_loss_mat_temp = matrix(0, idx_intcpt, n_resps)
    grad_loss_mat_temp[, idx_resps_conti] <- grad_loss_conti_mat
    grad_loss_mat_temp[, idx_resps_multi] <- grad_loss_multi_mat
    grad_loss_mat_temp[, idx_resps_ordin] <- grad_loss_ordin_mat
    grad_loss_mat = grad_loss_mat_temp * white_coefs
    
    pp = Se %*% (W1^2) %*% t(Sr)
    E = expm(pp)
    h = sum(diag(E)) - n_nodes
    SetEtSr = t(Se) %*% t(E) %*% Sr
    grad_h_mat_coef = SetEtSr * (2 * W1)
    grad_h_mat = matrix(0, idx_intcpt, n_resps)
    grad_h_mat[-idx_intcpt, ] = grad_h_mat_coef 
    grad_mat = grad_loss_mat + grad_h_mat * (rho1 * h + alpha)
    grad_vec  = vectorizing(grad_mat, idx_intcpt)
    (grad1_r = grad_vec + rho * (w - beta + mu))
    
    sourceCpp("./grad0_mix.cpp")
    
  }
  
  ## Wrapping grad function ####
  grad1 = function(w, alpha, rho1){
    grad0cpp_mix(w, beta, mu, rho, alpha, rho1, X, Y, white_coefs,
                 rev_interv_info, rev_interv_info_design, n_obs_by_node_conti,
                 idx_intcpt, n_resps, n_obs, n_nodes, n_coefs, 
                 as.numeric(idx_resps_conti), as.numeric(idx_resps_multi), 
                 as.numeric(idx_resps_ordin),
                 as.numeric(idx_multi_nodes), as.numeric(idx_ordin_nodes),
                 Sr, Se, SSr, as.numeric(rj_lev), as.numeric(rj_lev_sub), 
                 T_H, T_Mu, Y_tilde)}
  
  if (F){
    grad1_cpp = grad1(w, alpha, rho1)
    cbind(grad1_r, grad1_cpp)
  }
  
  
  while (rho1 < rho1_max) {
    
    if (verbose == 1){ cat("rho 1 =", rho1, " and alpha = ", alpha, "\n") }
    
    
    if (ver == "1"){    
      obj2 = function(w){obj1(w, alpha, rho1)}
      grad2 = function(w){grad1(w, alpha, rho1)}
      
    } else {
      sigma_sq_hat = (1/n_obs) * 
        apply(Y[,idx_resps_conti] - X%*%W[,idx_resps_conti], 2, function(x) sum(x^2))
      obj2 = function(w){obj1(w, sigma_sq_hat, alpha, rho1)}
      grad2 = function(w){grad1(w, sigma_sq_hat, alpha, rho1)}
      
    }
    
    
    if (opt_method == "lbfgs"){
      
      w_new_temp = lbfgs::lbfgs(obj2, grad2, w, invisible = 1)$par
      
    } else if (opt_method =="lbfgs-armijo") {
      
      w_new_temp = lbfgs::lbfgs(obj2, grad2, w, invisible = 1, 
                                linesearch_algorithm = 
                                  "LBFGS_LINESEARCH_BACKTRACKING_ARMIJO")$par
      
    } else if (opt_method == "nloptr") {
      
      w_new_temp = nloptr::lbfgs(w, obj2, grad2)$par
      
    }
    
    # W_new_temp = matricizing(w_new_temp, n_expls, n_resps)
    # W_new = W_new_temp * white_coefs
    # 
    # W1_new = W_new[-idx_intcpt, ] # coef part only
    # h_new = sum(diag(expm::expm(Se %*% W1_new^2 %*% t(Sr)))) - n_nodes

    h_new = arma_h_cal(w_new_temp, white_coefs, idx_intcpt, n_resps, n_coefs,
                       Sr, Se, n_nodes)
    
    if (verbose == 1){ cat("h_new =", h_new, "\n") }
    if (verbose == 1){ cat("convergence=", opt_result$convergence, "\n") }
    if (verbose == 1){ cat("=======================================\n") }
    
    if (h_new > 0.25 * h) rho1 = rho1 * 10
    
    # w = vectorizing(W_new, idx_intcpt) 
    w = w_new_temp
    h = h_new
    
    if (h <= h_tol) break
    
    alpha <- alpha + rho1 * h
    
  }
  
  W_new = matricizing(w, n_expls, n_resps)
  
  return(W_new)
  
}


