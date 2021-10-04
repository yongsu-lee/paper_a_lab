#######################################################
## Table 3 (Small Graphs, Nom vs Ord, Results) ####
#######################################################

if(T){ # for debugging, delete when submitting
  
  queue_arg = read.table("./codes/functions/queue_list_simu2", sep =",", 
                         strip.white = T)
  queue_args = queue_arg[1,]
  graph_type = as.character(queue_args[1])
  method = as.character(queue_args[2])
  iter = as.integer(queue_args[3])
  n_lams = 3 # the number of tuning parameters
  save_file_name = paste0("./results/simu2_small_results/simu2_small_",
                          graph_type,"_",method,"_",iter,".rds")
  
}

if(F){ # for server running
  
  # Clear the memory
  rm(list=ls())
  
  # Load required packages and functions
  source("00_load_ftns.R")
  
  # Set seed generating system
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  
  #!/usr/bin/env Rscript
  queue_args = commandArgs(trailingOnly=TRUE)
  graph_type = as.character(queue_args[1])
  method = as.character(queue_args[2])
  iter = as.integer(queue_args[3])
  save_file_name = paste0("simu2_small_",graph_type,"_",method,"_",iter,".rds")
  
  n_lams = 30
}


## Set up parameters
seed_para = 1
n_obs = 50 # the number of observations
n_conti = 6 # the number of continuous nodes
n_multi = 4 # the number of multinomial nodes
n_nodes = n_conti + n_multi # the number of total nodes


## Generate graphs, parameters and datapoints 

# Randomly assign continuous or multinomial nodes
types_by_node = gen_node_types(n_conti, n_multi, 0, seed = seed_para)

# Generate the number of levels for each node randomly 
n_levels_by_node = gen_node_levels(types_by_node, 4, seed = seed_para)
#.. it would be 1 if the node is continuous  

# Generate the true graph and obtain corresponding adjacency matrix
graph_set = gen_graph_adj_mat(n_nodes, graph_type, seed = seed_para)
graph_true = graph_set$graph_true
A_true = graph_set$A_true

# Generate the grand parameter matrix W 
W_true = gen_para(A_true, types_by_node, n_levels_by_node, seed = seed_para)

# Generate datapoints
data_input = gen_data(n_obs, A_true, graph_true, W_true, seed = iter)
#.. Note that the seed set with the 'iter' number


## Discretizing continuous nodes either 2 or 4 levels
conti_nodes = which(types_by_node == "c")
two_level_nodes = conti_nodes[1:(n_conti/2)]
four_level_nodes = conti_nodes[(n_conti/2+1):n_conti]

# Discretizing into 2 levels
temp = data_input[two_level_nodes]
temp1 = lapply(temp, function(x) 
  cut(x, c(min(x)-1, median(x), max(x)), labels = 1:2))

# Discretizing into 4 levels
temp = data_input[four_level_nodes]
temp2 = lapply(temp, function(x) 
  cut(x, c(min(x)-1, quantile(as.matrix(x), probs = c(0.25, 0.5, 0.75, 1))), 
      labels = 1:4))

discretized_data = data.frame(temp1, temp2)


## Generate a different type of dataset according to the method
data_input = 
  switch(method,
         # Oringinal data
         "mc" = { data_input },
         # Consider discretized nodes as 'nominal-level' nodes
         "mm" = { data_input[conti_nodes] <- discretized_data; data_input},
         # Consider discretized nodes as 'ordinal-level' nodes
         "mo" = {
           ordered = lapply(discretized_data, function(x) factor(x, ordered = T))
           data_input[conti_nodes] <- ordered; data_input
         } )


## Run main functions
result = glmdag(data_input, n_lams = n_lams, verbose = T,
                admm_args = admm_arg_ctrl(warm_start = F))


## Save the result
saveRDS(result, file = save_file_name)


## Collect the results

simu_title = "simu2_small"
subdir = "./results/simu2_small_cold/" # cold start case
# subdir = "./results/simu2_small/" # warm start case
methods =  c("mc","mm","mo")
graph_types = c("rand","sf","bi","sw","tree")
n_iter = 100

# thre=0.05  thre=0.1  thre=0.2  thre=0.3  thre=0.4  thre=0.5    thre=1 
thre_idx = 3
eval_crit = c("P","R", "E", "M","FP", "TPR", "FDR", "SHD", "JI")
thre_list = c(0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 1)

for (g in graph_types){ # g = "rand"
  
  graph_set = gen_graph_adj_mat(n_nodes, g, seed = seed_para)
  A_true = graph_set$A_true
  output_name = paste0(simu_title, "_", g, ".txt")
  write(x=NULL, file= output_name, append = F)
  
  for (m in methods){ # m = "mc"
    
    result_collect = array(dim = c(length(thre_list), length(eval_crit), n_iter))
    time_table = c()
    cat(paste0("Method = ", m), sep = "\n", append = T, file = output_name)
    
    for (i in 1:n_iter){ # i = 1
      
      file_name = paste0(subdir, simu_title, "_", g, "_", m, "_", i, ".rds")
      result = readRDS(file_name)
      
      A_est = result$pushed_A_est_by_lam
      sel_lam = result$sel_tun
      
      result_collect[,,i] <- as.matrix(
        eval_table(A_est, A_true, "obs", sel_lam, thre_list)[,-c(10,11)] )
      
    }
    
    ave_result = apply(result_collect, c(1,2), function(x) mean(x))
    colnames(ave_result) <- eval_crit
    ave_result_with_thre = cbind(thre = thre_list, ave_result)
    capture.output( print(ave_result_with_thre, print.gap = 3), 
                    append = T, file = output_name )
    # cat(paste0("Time: ", round(mean(time_table),2), "sec."), sep = "\n\n", 
    #     append = T,  file = output_name)
  }
  
}

