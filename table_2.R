#######################################################
## Table 2 (CD vs GLMDAG results) ####
#######################################################


if(T){ # for debugging, delete when submitting
  
  queue_arg = read.table("./codes/functions/queue_list_simu1", sep =",", 
                         strip.white = T)
  queue_args = queue_arg[1,]
  graph_type = as.character(queue_args[1])
  iter = as.integer(queue_args[2])
  n_lams = 3 # the number of tuning parameters
  save_file_name = paste0("./results/simu1_n_50_results/simu1_n_50_",
                          graph_type,"_",iter,".rds")
  
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
  iter = as.integer(queue_args[2])
  save_file_name = paste0("simu1_n_50_",graph_type,"_",iter,".rds")
  
}

## Generate graphs, parameters and datapoints 

# Set up parameters
seed_para = 1
n_obs = 50 # the number of observations
n_nodes = 30 # the number of total (multinomial) nodes

# Generate the true graph and obtain corresponding adjacency matrix
graph_set = gen_graph_adj_mat(n_nodes, graph_type, seed = seed_para)
graph_true = graph_set$graph_true
A_true = graph_set$A_true
n_edge = sum(A_true) 

# Set all the nodes as multinomial nodes
types_by_node = rep("m", n_nodes)

# Generate the number of levels for each node randomly 
n_levels_by_node = gen_node_levels(types_by_node, 4, seed = seed_para)
#.. it would be 1 if the node is continuous  

# Generate the grand parameter matrix W 
W_true = gen_para(A_true, types_by_node, n_levels_by_node, seed = seed_para)

# Generate datapoints
data_input = gen_data(n_obs, A_true, graph_true, W_true, seed = iter)
#.. Note that the seed set with the 'iter' number

# Convert dataset for cd algorithm
data_input_cd = conv_to_cd_data(data_input) 

## Generate result containers

## Run main functions - GLMDAG
result_glmdag = glmdag(data_input, n_lams = n_lams, verbose = T)

## Run main functions - CD
result_cd = cd.run(indata = data_input_cd)

result <- list(glmdag = result_glmdag, cd = result_cd)

## Save the result
saveRDS(result, file = save_file_name)