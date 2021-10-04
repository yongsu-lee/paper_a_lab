## Setting ####

rm(list=ls())
RNGkind("Mersenne-Twister", "Inversion", "Rejection")

if (Sys.info()[1] == "Darwin"){
  
  source("00_load_ftns.R")

  # queue_arg = read.table("./queue_list", sep =",", strip.white = T)
  # queue_args = queue_arg[1,]
  # graph_type = as.character(queue_args[1])
  # iter = as.integer(queue_args[2])
  
} else {
  
  source("00_load_ftns.R")
  
  #!/usr/bin/env Rscript
  queue_args = commandArgs(trailingOnly=TRUE)
  # graph_type = as.character(queue_args[1])
  # iter = as.integer(queue_args[2])
  
}

## Set up parameters ####
seed_para = 2
seed_iter = 1
n_obs = 20 # the number of observations
n_nodes = 10 # the number of total (multinomial) nodes
graph_type = "rand"
# types_by_node = rep("c", n_nodes)
types_by_node = gen_node_types(3, 3, 4, seed_para)
n_levels_by_node = gen_node_levels(types_by_node, 4, seed = seed_para)

## Generate a true graph
graph_set = gen_graph_adj_mat(n_nodes, graph_type, seed = seed_para)
graph_true = graph_set$graph_true
A_true = graph_set$A_true
n_edge = sum(A_true) 

## Generate a parameter matrix
(W_true = gen_para(A_true, types_by_node, n_levels_by_node, 
                   intcpt = "for_child_only", seed = seed_para))

## Generate a data matrix
(data_input = gen_data(n_obs, A_true, graph_true, W_true, seed = seed_iter))






