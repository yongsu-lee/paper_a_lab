gen_para = function(A_true, types_by_node = NULL, n_levels_by_node = NULL,
                    intcpt = c("always","for_child_only","none"),
                    seed = 1, range = c(0.5,1), prob_sparse = 0, 
                    ordin_pred_as = "num"){
  
  if (FALSE){
    ordin_pred_as = "num"
    intcpt = "for_child_only"
    range = c(0.5,1)
    prob_sparse = 0
    seed = 1
  }
  
  if (missing(intcpt)) {
    intcpt = "always"
  } else {
    intcpt = match.arg(intcpt)
  }
  
  graph_temp = igraph::graph_from_adjacency_matrix( A_true )
  
  if (!igraph::is_dag(graph_temp)) {
    stop("Provided adjacency matrix does not represent a DAG. Check the A_true.")
  }
  
  if (is.null(n_levels_by_node)){
    stop("n_levels_by_node is required!")
  }
  
  if (any(n_levels_by_node[(types_by_node == "m")] == 1)){
    stop("The number of levels for categorical nodes should be more than 1")
  }
  
  if (any(n_levels_by_node[(types_by_node == "o")] == 1)){
    stop("The number of levels for categorical nodes should be more than 1")
  }
  
  n_nodes = length(types_by_node)
  
  # n_resps_by_node = ifelse(types_by_node == "c", 
  #                          n_levels_by_node, n_levels_by_node-1)
  # .. Changed for symmetric lab ....
  n_resps_by_node = ifelse(types_by_node == "c", 
                           n_levels_by_node, n_levels_by_node)
  
  idx_resps_by_node = gen_indices_by_node(n_resps_by_node)
  
  # n_expls_by_node = n_resps_by_node
  # .. Changed for symmetric lab ....
  n_expls_by_node = ifelse(types_by_node == "c", 
                           n_levels_by_node, n_levels_by_node-1)
  
  if (ordin_pred_as == "num"){
    n_expls_by_node[types_by_node == "o"] = 1
  }
  
  idx_expls_by_node = gen_indices_by_node(n_expls_by_node)
  
  n_row = sum(n_expls_by_node)
  n_col = sum(n_resps_by_node)
  
  W_temp = matrix(0, n_row, n_col)
  
  set.seed(seed)
  for (j in 1:n_nodes){ #j = 1
    
    idx_resps_j = idx_resps_by_node[[j]]
    
    for (i in 1:n_nodes){ # i = 2
      
      idx_expls_i = idx_expls_by_node[[i]]
      
      if (A_true[i,j] == 1){
        
        n_cell = length(idx_expls_i) * length(idx_resps_j)
        W_temp[idx_expls_i, idx_resps_j] = sample(c(0, 1), n_cell, replace = T,
                                                  prob = c(prob_sparse, 1-prob_sparse))
      }
    }
  }
  
  # Adding intercept term ----
  if (intcpt == "for_child_only") {
    
    idx_intcpts = which(apply(W_temp, 2, sum) != 0)
    W_true = rbind(W_temp, rep(0,n_col))
    
    for (j in idx_intcpts) W_true[n_row + 1, j] = 1
    
  } else if (intcpt == "always") {
    
    W_true = rbind(W_temp, rep(1,n_col))
    
  } else {
    
    W_true = W_temp
  }
  
  # .. Changed for symmetric lab ....
  # Set the intercept of the first level of each node as zero (identifiability) ----
  first_level_intcpt = sapply(idx_resps_by_node, function(x) x[1])
  W_true[n_row + 1, first_level_intcpt] <- 0 
  
  # Assign generated coefficients  ----
  min_para = min(range); max_para = max(range)
  idx_nnz = which(W_true  == 1)
  n_nnz = sum(W_true)
  coefs = round(runif(n_nnz, min_para, max_para), 4)
  
  # Flipping the signs ----
  signs = sample(c(-1,1), n_nnz, replace = T)
  W_true[idx_nnz] = coefs * signs
  
  # .. Changed for symmetric lab ....
  # Put restriction across the level by level coefficients
  for (j in 1:n_nodes){ #j = 2
    
    idx_resps_j = idx_resps_by_node[[j]]
    
    for (i in 1:n_nodes){ # i = 6
      
      idx_expls_i = idx_expls_by_node[[i]]
      
      if (A_true[i,j] == 1){
        
        W_ij_transpose = t(W_true[idx_expls_i, idx_resps_j, drop = F])
        W_true[idx_expls_i, idx_resps_j] <- t(scale(W_ij_transpose))
        
      }
    }
  }
  
  
  # Attaching types of nodes to W_true ----
  attr(W_true, "types_by_node") = types_by_node
  attr(W_true, "n_levels_by_node") = n_levels_by_node
  if (any(types_by_node == "o")){
    attr(W_true, "ordin_pred_as") = ordin_pred_as 
  }
  attr(W_true, "n_resps_by_node") = n_resps_by_node
  
  colnames(W_true) <- resp_names(n_nodes, types_by_node, n_resps_by_node)
  
  if (intcpt == "none") {
    rownames(W_true) <- c(expl_names(n_nodes, types_by_node, n_expls_by_node,
                                     ordin_pred_as = ordin_pred_as))
    
  } else {
    rownames(W_true) <- c(expl_names(n_nodes, types_by_node, n_expls_by_node,
                                     ordin_pred_as = ordin_pred_as), "(Intcpt)")}
  
  
  return(W_true)
  
}
