rm(list = ls()); 

source("sourcing_nocpp.R")
sourceCpp("./arma_trexpm.cpp")
sourceCpp("./obj0_mix.cpp"); 
sourceCpp("./grad0_mix.cpp"); 
sourceCpp("./lkhd_mix.cpp")
sourceCpp("./lkhd_multi.cpp")


if(F){
  queue_arg = read.table("queue_list", sep =",", strip.white = T)
  queue_args = queue_arg[1,]
}

#!/usr/bin/env Rscript
queue_args = commandArgs(trailingOnly=TRUE)
graph_type = as.character(queue_args[1])
method = as.character(queue_args[2])
iter = as.integer(queue_args[3])

save_file_name = paste0("simu1_large_",graph_type,"_",method,"_",iter,".rds")

## generate graph and data ####

n_lams = 30
n_conti = 16
n_multi = 14
n_obs = 100

n_nodes = n_conti + n_multi 


seed_para = 1
graph_set = gen_graph_adj_mat(n_nodes, graph_type, seed = seed_para)
graph_true = graph_set$graph_true
A_true = graph_set$A_true

types_by_node = gen_node_types(n_conti, n_multi, 0, seed = seed_para)
n_levels_by_node = gen_node_levels(types_by_node, 4, seed = seed_para)

(W_true = gen_para(A_true, types_by_node, n_levels_by_node, seed = seed_para))
data_input = gen_data(n_obs, A_true, graph_true, W_true, seed = iter)

# discretizing according to the method

conti_nodes = which(types_by_node == "c")
two_level_nodes = conti_nodes[1:(n_conti/2)]
four_level_nodes = conti_nodes[(n_conti/2+1):n_conti]

temp = data_input[two_level_nodes]
temp1 = lapply(temp, function(x) cut(x, c(min(x)-1, median(x), max(x)), labels = 1:2))
temp = data_input[four_level_nodes]
temp2 = lapply(temp, function(x) 
               cut(x, c(min(x)-1, 
                       quantile(as.matrix(x), probs = c(0.25, 0.5, 0.75, 1))), 
                  labels = 1:4))

discretized_data = data.frame(temp1, temp2)

data_input = 
  switch(method,
         "mc" = { data_input },
         "mm" = { data_input[conti_nodes] <- discretized_data; data_input},
         "mo" = {
           ordered = lapply(discretized_data, function(x) factor(x, ordered = T))
           data_input[conti_nodes] <- ordered; data_input
           } )

tic <- proc.time()
result = noteargis(data_input, n_lams = n_lams, verbose = T)
toc <- proc.time() - tic

simu_result = list(result = result, timing = toc)

saveRDS(simu_result, file = save_file_name)
