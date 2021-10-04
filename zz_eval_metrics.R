##########################################################
## Script for generating Evaluation Related Functions ####
##########################################################

## Sub Functions and package

skel = function(A){
  result_temp = A + t(A)
  A_skeleton = (result_temp != 0) * 1
  return(A_skeleton)
}

rev_edge = function(A_true, A_est, type = "interv"){
  
  if (type == "interv"){
    
    rev_edges = sum((A_est == 1) & (t(A_true) == 1))
    return(rev_edges)
    
  }
  
  if (type == "obs"){
    
    cpdag_A = dag2cpdag(A_true)
    cpdag_est_A = dag2cpdag(A_est)
    
    # find out undirected edges from true/estiamte graphs
    comp1 = cpdag_A == t(cpdag_A)
    comp2 = cpdag_est_A == t(cpdag_est_A)
    
    rev_edges_temp = ((A_est==1) & (t(A_true)==1))
    
    # Among the rev_edges_temp, we will exclude case that
    # ... corresponding edge is undirected in both
    # ... CPDAGs of estimated and true graph
    non_R = sum((rev_edges_temp & comp1) & comp2)
    return(sum(rev_edges_temp) - non_R)
    
  }
}


## -----------------------------

## Main Functions

# (1) Save several criteria across the path ----

eval_by_lam = function(A_est_by_lam, A_true, n_edges, type = "interv"){
  
  P = R = E = M = FP = SHD = JI = c()
  
  n_lams = length(A_est_by_lam)
  
  for (ell in 1:n_lams){ # ell = 2
    
    A_est = A_est_by_lam[[ell]]
    
    P[ell] = sum(A_est)
    if (ell > 1 & P[ell] == 0){P[ell] <- 1e-8} 
    R[ell] = rev_edge(A_true, A_est, type)
    
    skel_A_true = skel(A_true)
    E[ell] = sum((skel_A_true == 1) & (A_est == 1)) - R[ell]
    
    if (ell > 1 & E[ell] == 0){E[ell] <- 1e-8} 
    
    M[ell] = n_edges - E[ell] - R[ell] 
    FP[ell] = P[ell] - R[ell] - E[ell]
    
    SHD[ell] = R[ell] + M[ell] + FP[ell]
    JI[ell] = E[ell] / (P[ell] + n_edges - E[ell])
    
  }
  
  func_result = list(P = P, E = E, JI = JI, R = R, M = M, FP = FP, SHD = SHD)
  return(func_result)
  
}

# (2) Select Model using the result from eval_by_lam ----

mod_sel = function(one_simu_result, n_edges, sel_crit = "JI"){
  
  sel_idx = ifelse (sel_crit == "JI",  which.max(one_simu_result$JI), 
                    which.min(one_simu_result$SHD))
  
  P= one_simu_result$P[sel_idx]
  E = one_simu_result$E[sel_idx]
  R = one_simu_result$R[sel_idx]
  M = one_simu_result$M[sel_idx]
  FP = one_simu_result$FP[sel_idx]
  TPR = E/n_edges
  FDP = (R + FP)/P
  SHD = one_simu_result$SHD[sel_idx]
  JI = one_simu_result$JI[sel_idx]

  func_result = data.frame(P, E, TPR, JI, R, M, FP, FDP, SHD)
  
  return(func_result)
  
}

eval_table = function(A_est_list, A_true, type, K_list, single = FALSE){
  
  if(F){ # debugging
    A_est_list = A_est; type = "obs"; K_list = sel_lam; 
    thre_list = c(0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 1)
    single = FALSE
  }

  thre_list = c(0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 1)
  
  if (single == FALSE) {
    n_thre = length(thre_list)
    n_edges = sum(A_true == 1)
  } else {
    n_thre = 1
    n_edges = sum(A_true == 1)
  }
  
  P = R = E = M = FP = SHD = JI = c()
  
  for (u in 1:n_thre){ # u = 1
    
    if (single == FALSE){
      K = K_list[u]
      A_est = A_est_list[[K]]
    } else {
      A_est = A_est_list
    }
    
    P[u] = sum(A_est)
    R[u] = rev_edge(A_true, A_est, type)
    
    skel_A_true = skel(A_true)
    E[u] = sum((skel_A_true == 1) & (A_est == 1)) - R[u]
    M[u] = n_edges - E[u] - R[u]
    FP[u] = P[u] - R[u] - E[u]
    
    TPR = E/n_edges
    FDP = (R + FP)/P
    SHD[u] = R[u] + M[u] + FP[u]
    JI[u] = E[u] / (P[u] + n_edges - E[u])
    
  }
  
  if (single == FALSE){
    eval_result = data.frame(P = P, E = E, TPR = TPR,  JI = JI, 
                             R = R, M = M, FP = FP, FDP = FDP, SHD = SHD, 
                             K = K_list, thre = thre_list)
  } else {
    eval_result = data.frame(P = P, E = E, TPR = TPR,  JI = JI, 
                             R = R, M = M, FP = FP, FDP = FDP, SHD = SHD)
  }
  
  return(eval_result)
  
}
