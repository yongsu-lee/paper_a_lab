gen_data = function(n_obs, A_true, graph_true, W_true, ordin_pred_as = "num",
                    interv_info=NULL, seed=NULL, df = TRUE) {
  
  if(FALSE){ # debugging input
    interv_info = matrix(FALSE, n_obs, ncol(A_true))
    # for( j in 0:(n_nodes-1)){interv_info[(10*j + 1):(10*(j+1)),(j+1)] <- TRUE }
    rev_interv_info = !interv_info
    interv_info = matrix(FALSE, n_obs, n_nodes)
    ordin_pred_as = "num"
    seed = 1; df = TRUE
  }
  
  if (missing(n_obs)) stop("Number of observation?")
  if (missing(A_true)) stop("True adjacency matrix?")
  if (missing(graph_true)) stop("True graph(igraph format)?")
  if (missing(W_true)) stop("True paramter?")
  
  n_nodes = ncol(A_true)
  
  if (is.null(interv_info)){
    interv_info = matrix(FALSE, n_obs, n_nodes)
  }
  
  nrow_interv_info = nrow(interv_info)
  ncol_interv_info = ncol(interv_info)
  
  if (n_obs != nrow_interv_info) stop("'n_obs' is not consistnet with the interv_info")
  if (n_nodes != ncol_interv_info) stop("Check the number of columns of interv_info")
  
  types_by_node = attr(W_true, "types_by_node")
  n_levels_by_node = attr(W_true, "n_levels_by_node")
  
  # n_resps_by_node = ifelse(types_by_node == "c", 
  #                          n_levels_by_node, n_levels_by_node-1)
  # .. Changed for symmetric lab ....
  n_resps_by_node = n_levels_by_node
  
  # n_expls_by_node = n_resps_by_node
  # .. Changed for symmetric lab ....
  n_expls_by_node = ifelse(types_by_node == "c", 
                           n_levels_by_node, n_levels_by_node-1)
  
  if (ordin_pred_as == "num"){
    n_expls_by_node[types_by_node == "o"] = 1
  }
  
  idx_resps_by_node = gen_indices_by_node(n_resps_by_node)
  idx_expls_by_node = gen_indices_by_node(n_expls_by_node)
  
  parents_by_node = lapply(as.data.frame(A_true), function(x) which(x==1))
  n_obs_by_node = rep(0, n_nodes)
  
  if (any(rownames(W_true) == "(Intcpt)")){
    n_expls = nrow(W_true) - 1
    idx_intcpt = nrow(W_true)
  } else {
    n_expls = nrow(W_true)
    idx_intcpt = NULL
  }
  
  X1_temp = matrix(NA, n_obs, n_expls)
  Z = matrix(NA, n_obs, n_nodes)
  
  topo_info = topo_sort(graph_true)
  
  set.seed(seed)
  # seq_len(n_nodes)
  for (j in seq_len(n_nodes)){ # j = 6
    
    resp_node = topo_info[j]
    q_j = n_levels_by_node[resp_node]
    node_type = types_by_node[resp_node]
    idx_resps_j = idx_resps_by_node[[resp_node]]
    n_resps_j = length(idx_resps_j)
    
    expl_nodes = parents_by_node[[resp_node]]
    
    ############################
    ## if it is a root node #### 
    ############################
    
    ## Note: Intervention processing is pointless for root nodes.
    
    if (length(expl_nodes) == 0L){ 
      
      ####################
      ## categorical cases (nominal/ordinal)
      
      if (node_type == "m" | node_type == "o"){
        
        Z_j = try_gen_levels(q_j, n_obs)
        
        ## .. Intervention processing is pointless for root nodes.
        interv_j = interv_info[ ,resp_node, drop = F]
        n_obs_interv_j = sum(interv_j)
        n_obs_by_node[resp_node] = n_obs - n_obs_interv_j
        
        Z[, resp_node] = Z_j
        X1_temp[ , idx_expls_by_node[[resp_node]]] <- 
          gen_X1_temp(node_type, Z_j, q_j, ordin_pred_as)
        
        n_obs_by_node[resp_node] = n_obs - n_obs_interv_j
        
        ##########################################
        ## Gaussian case (or continuous, possibly)
        
      } else {
        
        Z_j = matrix(rnorm(n_obs * n_resps_j), ncol = n_resps_j )
        
        ## .. Intervention processing is pointless for root nodes.
        interv_j = interv_info[ ,resp_node, drop = F]
        n_obs_interv_j = sum(interv_j)
        n_obs_by_node[resp_node] = n_obs - n_obs_interv_j
        # Z_j[interv_j,] <- rnorm(n_obs_interv_j)
        
        Z[, resp_node] = Z_j
        X1_temp[, idx_expls_by_node[[resp_node]]] = Z_j
        
      }
      
      ###########################################################
      ## if it is not a root node (at least, has one parent) #### 
      ###########################################################
      
    } else { 
      
      idx_expls_j = unlist(idx_expls_by_node[expl_nodes])
      n_expls_j = length(idx_expls_j)
      
      ##########################################
      ## categorical cases 
      if (node_type == "m" | node_type == "o"){
        
        etas1 = X1_temp[, idx_expls_j, drop = F] %*% 
          (W_true[idx_expls_j, idx_resps_j, drop = F])
        
        ## Add intercepts
        if (any(rownames(W_true) == "(Intcpt)")){
          etas = sweep(etas1, 2, W_true[idx_intcpt, idx_resps_j, drop = F], "+")
        } else {
          etas = etas1
        }
        
        ## Calculate Probs under baseline-category logit
        Mu = switch( node_type, 
                     "m" = {
                       exp(etas) },
                     "o" = {
                       H = etas %*% lower.tri(matrix(1,q_j,q_j), diag=T)
                       exp(H) } )
        
        # temp_deno = 1 + cbind(rowSums(Mu))
        # probs = cbind(sweep(Mu, 1, temp_deno, "/"), 1/temp_deno)
        # .. Changed for symmetric lab ....
        probs = sweep(Mu, 1, cbind(rowSums(Mu)) , "/")
        
        ## Intervention data
        interv_j = interv_info[ , resp_node, drop = F]
        n_obs_interv_j = sum(interv_j)
        n_obs_by_node[resp_node] = n_obs - n_obs_interv_j
        probs[interv_j, ] <- 1/q_j 
        # .. Do not use 'drop=F' option for integrating NULL indexing.
        
        Z_j = try_gen_levels(q_j, n_obs, probs)
        Z[, resp_node] = Z_j
        X1_temp[, idx_expls_by_node[[resp_node]]] <- 
          gen_X1_temp(node_type, Z_j, q_j, ordin_pred_as)
        
        ## Gaussian case (or continouos, possibly)            
      } else { 
        
        Z_j = X1_temp[, idx_expls_j, drop = F] %*% 
          (W_true[idx_expls_j, idx_resps_j, drop = F]) + 
          matrix(rnorm(n_obs * n_resps_j ), ncol = n_resps_j )
        
        ## intervention data
        interv_j = interv_info[,resp_node, drop = F]
        n_obs_interv_j = sum(interv_j)
        n_obs_by_node[resp_node] = n_obs - n_obs_interv_j
        # Z_j[interv_j,] <- 10
        Z_j[interv_j,] <- rnorm(n_obs_interv_j)
        
        X1_temp[, idx_expls_by_node[[resp_node]]] = Z_j
        Z[, resp_node] = Z_j
        
      }
      
    } # end of if-else sttt of length(expl_nodes) == 0L
    
  } # end of for loop (j in 1:n_nodes)
  
  if (df == T){
    Z <- as.data.frame(Z)
    Z[, (types_by_node=="m") | (types_by_node=="o") ] <- 
      lapply(Z[,types_by_node=="m" | types_by_node=="o" , drop = F], as.factor)
  }
  
  colnames(Z) <- paste0("Z", 1:n_nodes)
  attr(Z,"n_obs_by_node") = n_obs_by_node
  
  return(Z)
  
}


try_gen_levels = function(q_j, n_obs, probs = NULL, trials = 100){
  
  k = 0;
  
  while (k < trials){
    
    if (is.null(probs)){
      Z_j = factor(sample(0:(q_j-1), n_obs, replace = T))
    } else {
      Z_j = factor( apply(probs, 1, function(x) 
        sample(0:(q_j-1), 1, replace=T, prob = x)) )
    }
    
    if (q_j == nlevels(Z_j)) break
    
    k = k + 1
    if (k > trials){
      stop("Something's wrong. Single level column has been generated.")}
  } 
  
  return(Z_j)
  
}

gen_X1_temp = function(node_type, Z_j, q_j = NULL, ordin_pred_as = NULL){
  
  if (is.null(ordin_pred_as)) ordin_pred_as <- "num"
  if (ordin_pred_as == "nom") node_type <-  "m"
  
  switch(node_type,
         "c" = {
           X1_temp_j <- Z_j },
         "m" = {
           X1_temp_j <-  
             model.matrix(~factor(Z_j) + 0)[, -q_j, drop = F] },
         "o" = {
           X1_temp_j <-  as.numeric(Z_j) }
  )
  
  return(as.matrix(X1_temp_j))
}

### Check Up List ####

# - deal with intercept part (just in case it does not exist)
# - double check for gaussian case (especially for intervention part)