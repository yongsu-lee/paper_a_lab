#' Control ADMM Arguments
#'
#' @param max_iter a positive integer, the number of loops in ADMM algorithm. Default is \code{100}.
#' @param rho a positive double scalar, a penalized parameter in the ADMM loop. Default is \code{1}.
#' @param kappa a postive double scalar, a relaxation pameter in the ADMM loop for the better convergnece. Recommended value is between 1.6 and 1.8(Default).
#' @param abs_tol a positive double scale, an absolute tolerance for stopping criterion in the ADMM loop.
#' @param rel_tol positive double scale, an relative tolerance for stopping criterion in the ADMM loop.
#' @param eps a positive scalar, a tolerance of estimates of parameter for the various purpose. Default is \code{1e-8}
#' @param 
#'
#' @return a list object with arguments for customizing ADMM algorithm.
#' @export
#'
admm_arg_ctrl <- function(max_iter = 100, rho = 1.0, kappa = 1.8, abs_tol = 1e-4, 
                          rel_tol = 1e-2, eps = 1e-8, warm_start = T, 
                          inner_verbose = F){
  
  value = list(max_iter = max_iter, rho = rho, kappa = kappa, abs_tol = abs_tol,
               rel_tol = rel_tol, eps = eps, warm_start = warm_start, 
               inner_verbose = inner_verbose)
  
  invisible(value)
  
}

layout_set = function(graph_true, graph_type){
  
  if (graph_type == "bi") {
    out = layout_as_bipartite(graph_true)
  } else if (graph_type == "sw") {
    out = layout_in_circle(graph_true)
  } else if (graph_type == "tree"){
    out = layout_as_tree(graph_true)
  } else if (graph_type == "sf"){
    out = layout.fruchterman.reingold(graph_true)
  } else if (graph_type == "rand"){
    out = layout_nicely(graph_true)
  } else {
  }
  
  return(out)
}

fancy_dag_plot = function(A_true, types_by_node, black_list=c(),
                          font_size = 1, title = "DAG graph", cex.main = 2, 
                          legend = T, node_names = "num", ...){
  
  if (node_names == "num"){
    num_names = 1:ncol(A_true)
  }
  
  graph_true = igraph::graph_from_adjacency_matrix(A_true)
  graph_true$name <- paste0("True Graph (graph_type = ", attr(graph_true, "graph_type"),")")
  V(graph_true)$color<- ifelse(types_by_node == "m", "black", "white")
  n_nodes = length(types_by_node)
  
  plot(graph_true,
       vertex.label = num_names,
       vertex.shape = ifelse(1:5 %in% black_list, "square", "circle"),
       vertex.label.color = ifelse(types_by_node == "m", "white", "black"),
       vertex.label.cex = font_size, 
       vertex.label.font = 2, 
       ...)
  title(main = title, cex.main = cex.main)
  if (legend == T){
  legend("topleft", legend=c("Multi-level", "Gaussian"),
         pch = c(19,1) , cex=1)
  }
}

gen_node_types = function(n_conti_nodes = 0, n_multi_nodes= 0,  n_ordin_nodes = 0,
                          seed = 1){
  if(F){ # debugging input
    n_conti_nodes=2; n_multi_nodes=3; n_ordin_nodes=1; seed = 1
  } # Done/Checked : 02/01/21  
  
  temp = c(rep("c", n_conti_nodes), rep("m", n_multi_nodes), rep("o", n_ordin_nodes))
  n_nodes = length(temp)
  
  set.seed(seed)
  types_by_node = sample(temp, n_nodes)
  
  names(types_by_node) = paste0("V", 1:n_nodes)
  return(types_by_node)
  
}

gen_node_levels = function(types_by_node, max_n_levels = 2, seed = 1){
  if(F){ # debugging input
    types_by_node = c("m","o","m","o","c")
    max_n_levels = 4
    seed = 1
  } # Done/Checked : 02/01/21   
  
  n_nodes = length(types_by_node)
  n_levels_by_node = rep(1, n_nodes)
  idx_categ_nodes = (types_by_node == "m" | types_by_node == "o")
  
  set.seed(seed)
  n_levels = sample(x = 1:(max_n_levels-1), 
                    size = sum(idx_categ_nodes), replace = TRUE) + 1
  
  n_levels_by_node[idx_categ_nodes] <- n_levels
  names(n_levels_by_node) <- paste0("V",1:n_nodes)
  
  return(n_levels_by_node)  
}

## For messages on algorithm iterations
## Called     : admm_loop(); *.W_update() 
## Last Updated : ??
ord = function(x) {
  
  ifelse(x == 1, "st", ifelse(x==2, "nd","th"))
  
} ##############################################################################


# Called        : 
# Last Updated  : 02/01/21
gen_indices_by_node = function(sizes_by_node){
  idx = split(1:sum(sizes_by_node),
              rep(seq_along(sizes_by_node), times = sizes_by_node))
  
  return(idx)
  
} ##############################################################################


## Called        : gen_para.R  
## Last Updated  : 02/06/21
resp_names = function(n_nodes, types_by_node, n_resps_by_node){
  
  
  letter = "Y"
  var_name_temp = rep(paste0(letter,1:n_nodes), times= n_resps_by_node)
  elmts_temp = as.list(n_resps_by_node)
  
  dummys_temp = as.list(rep(0,n_nodes))
  
  for (j in 1:n_nodes){ # j = 1
    
    node_type = types_by_node[j]
    node_n_resps = n_resps_by_node[j]
    
    if (node_type  == "m" | node_type == "o"){
      dummys_temp[[j]] <- paste0(rep(node_type, node_n_resps), ":D", 1:node_n_resps)
    } else if (node_type == "c") {
      dummys_temp[[j]] <- paste0(rep(node_type, node_n_resps), ":", 1:node_n_resps)
    }
  }
  
  dummys_temp[types_by_node=="c" & n_resps_by_node == 1] <- "c"
  dummys = unlist(dummys_temp)
  
  var_name = paste0(var_name_temp, dummys)
  
  return(var_name)
} ##############################################################################


## Called        : gen_para.R, shared_fun.R - gen_X1_design()
## Last Updated  : 02/06/21
expl_names = function(n_nodes, types_by_node, n_expls_by_node, ordin_pred_as){
  
  letter = "X"
  var_name_temp = rep(paste0(letter,1:n_nodes), times= n_expls_by_node)
  elmts_temp = as.list(n_expls_by_node)
  
  dummys_temp = as.list(rep(0,n_nodes))
  
  for (j in 1:n_nodes){ # j = 1
    
    node_type = types_by_node[j]
    node_n_expls = n_expls_by_node[j]
    
    if (node_type == "m" | node_type == "o"){
      dummys_temp[[j]] <- paste0(rep(node_type, node_n_expls), ":D", 1:node_n_expls)
    } else if (node_type == "c") {
      dummys_temp[[j]] <- paste0(rep(node_type, node_n_expls), ":", 1:node_n_expls)
    }
  }
  
  dummys_temp[types_by_node == "c" & n_expls_by_node == 1] <- ""
  dummys_temp[types_by_node == "o" & ordin_pred_as == "num"] <- ""
  dummys = unlist(dummys_temp)
  
  var_name = paste0(var_name_temp, dummys)
  
  return(var_name)
} ##############################################################################

## Called        : ordin.noteargis.R
## Last Updated  : 02/08/21
gen_n_xxxxx_by_node = function(role, types_by_node, n_levels_by_node, ordin_pred_as = "num"){
  
  n_resps_by_node = ifelse(types_by_node == "c", 
                           n_levels_by_node, n_levels_by_node-1)
  
  n_expls_by_node = n_resps_by_node
  if (ordin_pred_as == "num"){
    n_expls_by_node[types_by_node == "o"] = 1
  }
  
  if (role == "resps") { 
    return(n_resps_by_node) 
  } else if (role == "expls") {
    return(n_expls_by_node) 
  } else { 
    stop ("Incorrect 'role' input. Either 'resps' or 'expls.'")
  }
  
}

## Generate design matrix for data as predictors ###############################
## Called        : ordin.noteargis.R, noteargis.R
## Last Updated  : 03/31/21
gen_X1_design= function(data_input, types_by_node, n_expls_by_node, ordin_pred_as = "num"){
  
  if(F){
    types_by_node = rep("m",5)
    n_expls_by_node = gen_n_xxxxx_by_node("expls", types_by_node, n_levels_by_node)
    ordin_pred_as = "nom"
  }
  
  n_nodes = length(types_by_node)
  n_levels_by_node = sapply(data_input, nlevels)
  idx_expls_by_node = gen_indices_by_node(n_expls_by_node)
  
  temp_types_by_node = types_by_node
  if (ordin_pred_as == "num"){
    temp_types_by_node[temp_types_by_node == "o"] <- "c"
  } else if (ordin_pred_as == "nom") {
    temp_types_by_node[temp_types_by_node == "o"] <- "m"
  } else {
    stop ("Incorrect 'ordin_pred_as' input. Either 'nom' or 'num.'")
  }
  
  X1 = matrix(NA, dim(data_input)[1], 0)
  
  for (j in 1:n_nodes){ # j = 1
    
    q_j = n_levels_by_node[j]
    node_type = temp_types_by_node[j]
    
    x_temp = switch(node_type,
                    "c" = {
                      if (node_type == "c" & n_expls_by_node[j] > 1) {
                        data_input[, idx_expls_by_node[[j]]]
                      } else {
                        x_temp = data_input[,j] } },
                    "m" = {
                      x_temp = data_input[,j]
                      # model.matrix(~x_temp+0)[,-q_j]}
                      diag(q_j)[as.numeric(x_temp), -q_j] })
    
    X1 = cbind(X1, x_temp)
    
  }
  
  colnames(X1) <- expl_names(n_nodes, types_by_node, n_expls_by_node,
                             ordin_pred_as = ordin_pred_as)
  rownames(X1) <- 1:dim(data_input)[1]
  
  return(X1)
} ##############################################################################


## Generate design matrix for data as responses ################################
## Called        : ordin.noteargis.R, noteargis.R
## Last Updated  : 02/012/21
gen_Y_design= function(data_input, types_by_node, n_resps_by_node){
  
  n_nodes = length(types_by_node)
  idx_resps_by_node = gen_indices_by_node(n_resps_by_node)
  
  temp_types_by_node = types_by_node
  temp_types_by_node[types_by_node == "o"] <- "m"

  
  Y = matrix(NA, dim(data_input)[1], 0)
  
  for (j in 1:n_nodes){ # j = 1
    
    node_type = temp_types_by_node[j]
    
    y_temp = switch(node_type,
                    "c" = {
                      if (node_type == "c" & n_resps_by_node[j] > 1) {
                        data_input[, idx_resps_by_node[[j]]]
                      } else {
                        data_input[,j] } },
                    "m" = {
                      Z_j = data_input[,j]
                      q_j = nlevels(Z_j)
                      temp = diag(q_j)[as.numeric(Z_j), -q_j]}
    )
    
    Y = cbind(Y, y_temp)
  }
  
  colnames(Y) <- resp_names(n_nodes, types_by_node, n_resps_by_node)
  rownames(Y) <- 1:dim(data_input)[1]
  
  return(Y)
  
} ##############################################################################


## Check the input dataset 
## Called       : *.noteargis.R
## Last Updated : 02/08/21
## !! NOTE !! : this function should be updated for global purpose. Now, this ... 
### ... is only for the categorical.noteargis. 
check_input = function(data_input, types_by_node = NULL){
  
  if (is.data.frame(data_input)) {
    if (!all(sapply(data_input, class) == "factor")){
      idx_conti_nodes = (sapply(data_input, class) != "factor")
      message("There exist continuous variables. 
              These will be coerced to cateogrical(ordinal) data.")
      data_temp <- data_input
      data_input[, idx_conti_nodes] <- 
        as.data.frame(lapply(data_temp[,idx_conti_nodes], as.factor))
    } # otherwise, good to go. 
  } else {
    message("Input data is not a dataframe. It will be converted to a dataframe.")
    data_temp <- as.data.frame(data_input)
    data_input <- as.data.frame(lapply(data_temp, as.factor))
  }
  
  return(data_input)
} ##############################################################################



vectorizing = function(V, idx_intcpt){
  intcpts = cbind(V[idx_intcpt,])
  coefs = matrix(V[-idx_intcpt,], ncol = 1)
  v = rbind(coefs, intcpts)
  return(v)
}

matricizing = function(w, n_expls, n_resps, coef_only = F){
  n_coefs = n_expls * n_resps
  coefs = c(head(w, n = n_coefs))
  intcpts = c(tail(w, n = n_resps))
  
  W1 = matrix(coefs, nrow = n_expls, ncol = n_resps)
  
  if (coef_only == T) return(W1) else {
    W = rbind(W1, intcpts)
    rownames(W) <- NULL
    return(W)
  }
}

norm2 = function(x) {norm(as.matrix(x), "f")}

grp_shrinkage = function(x, a){
  z = max(0, 1- a/norm2(x)) * x
  return(z)
}

conv_to_cd_data = function(data){
  func_result = sparsebnUtils::sparsebnData(data, type = "discrete")
  return(func_result)
}

calc_A_est = function(Beta_new, data_info, eps = 1e-8){
  
  idx_intcpt = data_info$n_expls + 1
  Se = data_info$Se
  Sr = data_info$Sr
  
  Beta1_new = Beta_new[-idx_intcpt, ]
  
  K = Se %*% Beta1_new^2 %*% t(Sr)
  A_est = (K > eps)*1
  
  return(A_est)
}

calc_A_est_by_lam_cd = function(result_cd){
  
  n_lams = length(get.lambdas(result_cd))
  func_result = list(0)
  
  for (ell in 1:n_lams){ # ell = 1
    
    func_result[[ell]] = as.matrix(get.adjacency.matrix(result_cd[[ell]]))
    
  }
  
  return(func_result)
  
}

#' Print a \code{glmdag} object
#'
#' @param x fitted \code{glmdag} object
#' @param \dots additional arguments if available
#'
#' @return Display \code{graph_est_by_lam} and \code{A_est_by_lam}
#' @method print glmdag
#' @export
print.glmdag = function(x, ...){
  out = x
  cat("$graph_est_by_lam")
  print(out$graph_est_by_lam)
  cat("$A_est_by_lam")
  print(out$pushed_A_est_by_lam)
  cat("$Beta_est_by_lam")
  print(out$pushed_Beta_est_by_lam)
  
}

# print.glmdag_para = function(x, ...){
#   out = x
#   print(x, digits = 4)
# }

# print.glmdag_simu_data = function(x, ...){
#   out = x
#   options(digits = 2)
#   print(x)
#   options(digits = 7)
# }


## Sub Functions and package

skel = function(A){
  result_temp = A + t(A)
  A_skeleton = (result_temp != 0) * 1
  return(A_skeleton)
}

rev_edge = function(A_true, A_est, type = "obs"){
  
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

eval_by_lam = function(A_est_by_lam, A_true, n_edges, type = "obs"){
  
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


################################################################################
## (Soon be deprecated functions ) #############################################
################################################################################

## This function will be deprecated soon.
Y_names = function(n_nodes, types_by_node, n_resps_by_node){
  
  var_name_temp = rep(paste0("Y",1:n_nodes), times= n_resps_by_node)
  elmts_temp = as.list(n_resps_by_node)
  
  levels_temp = as.list(rep(0,n_nodes))
  
  for (j in 1:n_nodes){
    
    if (types_by_node[j] == "m"){
      levels_temp[[j]] <- paste0("=", 1:n_resps_by_node[j])
    } else if (types_by_node[j] == "c") {
      levels_temp[[j]] <- paste0(":", 1:n_resps_by_node[j])
    }
  }
  
  levels_temp[types_by_node=="c" & n_resps_by_node==1] <- ""
  levels = unlist(levels_temp)
  
  var_name = paste0(var_name_temp, levels)
  
  return(var_name)
  
} ##############################################################################



################################################################################
## (Old Eval Functions ) #######################################################
################################################################################

# skel = function(A){
#   result_temp = A + t(A)
#   A_skeleton = (result_temp != 0) * 1
#   return(A_skeleton)
# }
# 
# 
# rev_edge = function(A_true, A_est){
#   
#   cpdag_A = pcalg::dag2cpdag(A_true)*1
#   cpdag_est_A = pcalg::dag2cpdag(A_est)*1
#   
#   # find out undirected edges from true/estiamte graphs
#   comp1 = cpdag_A == t(cpdag_A)
#   comp2 = cpdag_est_A == t(cpdag_est_A)
#   
#   rev_edges_temp = ((A_est==1) & (t(A_true)==1))
#   
#   # Among the rev_edges_temp, we will exclude case that
#   # ... corresponding edge is undirected in both
#   # ... CPDAGs of estimated and true graph
#   non_R = sum((rev_edges_temp & comp1) & comp2)
#   return(sum(rev_edges_temp) - non_R)
# }

################################################################################