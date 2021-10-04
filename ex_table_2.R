#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
## Preamble ######################
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
rm(list=ls())

# install bioconductor packages if it not available in your system
if (!requireNamespace("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")}
BiocManager::install(c("graph", "RBGL"))

# install and load necessary packages
if (!require(devtools)){install.packages("devtools"); library(devtools)}
if (!require(igraph)){install.packages("igraph"); library(igraph)}
if (!require(discretecdAlgorithm)){install.packages("discretecdAlgorithm");
  library(discretecdAlgorithm)}

# load necessary packages
install_github("yongsu-lee/noteargis", force =T)
library(noteargis) # suggested methodology package

# install a dummy package to reproduce the simulation part.
# .. it contains all the simulation files, so it might takes a while.
install_github("yongsu-lee/noteargisEval", force =T)
library(noteargisEval)

# set seed generating system
RNGkind(kind="Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
## Section 5.4.1 Justification of Group-Structure based Algorithm ####
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# _ (1) Generate the Figure 1 ####

m = c(1,1,2,2,3,3,1,1,2,2,3,3,0,4,4,5,5,0,0,4,4,5,5,0)
M = matrix(m, nrow =4, ncol = 6, byrow = T); layout(mat = M)

n_nodes = 6
bi = attr(gen_adj_mat(n_nodes, "bi"), "graph_true")
rand = attr(gen_adj_mat(n_nodes, "rand"), "graph_true")
sf = attr(gen_adj_mat(n_nodes, "sf"), "graph_true")
sw = attr(gen_adj_mat(n_nodes, "sw"), "graph_true")
tree = attr(gen_adj_mat(n_nodes, "tree"), "graph_true")

pdf("figure1.pdf", width = 10)
layout(mat = M); par(mar=c(0,0.5,3,0.5)); par(oma=c(0,0,0,0))

dag_plot(bi, "c", 2, main = "Bipartite", layout = layout_as_bipartite(bi),
         vertex.label.font = 2, vertex.size = 40, edge.arrow.size = 0.5)

# dag_plot(rand, "c", 2, main = "Random DAG",
#          vertex.label.font = 2, vertex.size = 40, edge.arrow.size = 0.5)
## To generate the 'exactly' same graph in the paper, run below
layout_rand = readRDS(system.file("extdata", "layout_rand_table2.rds",
                                  package = "noteargis"))
dag_plot(rand, "c", 2, main = "Random DAG", layout = layout_rand,
         vertex.label.font = 2, vertex.size = 40, edge.arrow.size = 0.5)

# dag_plot(sf, "c", 2, main = "Scale-Free",
#          vertex.label.font = 2, vertex.size = 40, edge.arrow.size = 0.5)
## To generate the 'exactly' same graph in the paper, run below
layout_sf = readRDS(system.file("extdata", "layout_sf_table2.rds",
                                package = "noteargis"))
dag_plot(sf, "c", 2, main = "Scale-Free", layout = layout_sf,
         vertex.label.font = 2, vertex.size = 40, edge.arrow.size = 0.5)

dag_plot(sw, "c", 2, main = "Small-World", layout = layout_in_circle(sw),
         vertex.label.font = 2, vertex.size = 40, edge.arrow.size = 0.5)

dag_plot(tree, "c", 2, main = "Tree", layout = layout_as_tree(tree),
         vertex.label.font = 2, vertex.size = 40, edge.arrow.size = 0.5)
dev.off()


# _ (2) Generate the Table 2 ####

use_saved_data = 1
# if 1, it will used the saved results to generate table.
# if 0, it will run the functions (it will take much longer)

n_obs = 50
n_nodes = 6
types_by_node = rep("c", n_nodes)

iterations = 50
seed_para = 1
graph_types = c("bi","rand", "sf", "sw","tree")
methods = c("grp", "ungrp")
n_types = length(graph_types)
n_methods = length(methods)

n_resps_by_node = c(2,3,4,3,5,3)

total_grp_result_by_type = as.list(rep(0,n_types))
names(total_grp_result_by_type) <- graph_types
total_ungrp_result_by_type <- total_grp_result_by_type

for (graph_type in graph_types){ # graph_type = "bi"

  A_true = gen_adj_mat(n_nodes, graph_type, seed = seed_para)
  W_true = gen_para(A_true, types_by_node, n_resps_by_node, seed = seed_para)
  result_ungrp = c()
  result_grp = c()

  for (seed_data in 1:iterations){ # seed_data = 1

    if (use_saved_data == 1){
      obj_name = paste0("grp_", graph_type, "_data_",seed_data,"_table2.rds")
      fit.noteargis = readRDS(system.file("extdata", obj_name,
                                          package = "noteargisEval"))
    } else {

      data_input = gen_data(n_obs, A_true, W_true, seed = seed_data)
      fit.noteargis = noteargis(data_input, "c", n_resps_by_node, verbose = T)
    }

    Beta_est_by_lam_grp = fit.noteargis$Beta_est_by_lam
    result_grp = c(result_grp, min(sapply(Beta_est_by_lam_grp,
                                          function(x) norm(x-W_true,"f"))))

    if (use_saved_data == 1){
      obj_name = paste0("ungrp_", graph_type, "_data_", seed_data,"_table2.rds")
      fit.notears = readRDS(system.file("extdata", obj_name,
                                        package = "noteargisEval"))
    } else {
      fit.notears = noteargis(data_input, "c", verbose = T)
    }

    Beta_est_by_lam_ungrp = fit.notears$Beta_est_by_lam
    result_ungrp = c(result_ungrp, min(sapply(Beta_est_by_lam_ungrp,
                                              function(x) norm(x-W_true,"f"))))

  }

  total_ungrp_result_by_type[[graph_type]] <- result_ungrp
  total_grp_result_by_type[[graph_type]] <- result_grp

}

result = matrix(NA, ncol = n_types, nrow = 2)
colnames(result) <- graph_types; rownames(result) <- c("NOTEARGIS", "NOTEARS")


for (g in 1:n_types){ # g = 1
  graph_g = graph_types[g]
  result[1,g] <- mean(total_grp_result_by_type[[graph_g]])
  result[2,g] <- mean(total_ungrp_result_by_type[[graph_g]])
}

sink("table2.txt")
round(result,4)
sink()


#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
## Section 5.4.2 DAG Learning with Multinomial Nodes ####
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# _ (1) Generate the Figure 2 ####





# _ (2) Generate the Table 3 ####

use_saved_data = 1
# if 1, it will used the saved results to generate table.
# if 0, it will run the functions (it will take much longer)

n_obs = 50
n_nodes = 30
types_by_node = rep("m", n_nodes)

iterations = 50
seed_para = 1
graph_types = c("bi","rand", "sf", "sw","tree")
methods = c("noteargis", "cd")
n_types = length(graph_types)
n_methods = length(methods)

n_resps_by_node = gen_node_levels(types_by_node, max_n_levels = 4, seed = seed_para)

total_result_cd = as.list(rep(0,n_types))
names(total_result_cd) <- graph_types
total_result_noteargis <- total_result_cd
total_result_n_edges <- total_result_cd

sel_crit = "JI"

for (graph_type in graph_types){ # graph_type = "bi"

  A_true = gen_adj_mat(n_nodes, graph_type, seed = seed_para)
  n_edges = sum(A_true)

  W_true = gen_para(A_true, types_by_node, n_resps_by_node, seed = seed_para)
  result_cd = data.frame()
  result_noteargis = data.frame()

  for (seed_data in 1:iterations){ # seed_data = 1

    ## NOTEARGIS
    if (use_saved_data == 1){
      obj_name = paste0("noteargis_", graph_type, "_data_",seed_data,"_table3.rds")
      fit.noteargis = readRDS(system.file("extdata", obj_name,
                                          package = "noteargisEval"))
    } else {
      data_input = gen_data(n_obs, A_true, W_true, seed = seed_data)
      fit.noteargis = noteargis(data_input, "m", n_resps_by_node, verbose = T)
    }

    temp_noteargis_result = eval_by_lam(fit.noteargis$A_est_by_lam, A_true, n_edges)
    temp_noteargis_sel = mod_sel(temp_noteargis_result, n_edges, sel_crit = sel_crit)
    result_noteargis = rbind(result_noteargis, temp_noteargis_sel)

    ## CD
    if (use_saved_data == 1){
      obj_name = paste0("cd_", graph_type, "_data_",seed_data,"_table3.rds")
      fit.cd = readRDS(system.file("extdata", obj_name, package = "noteargisEval"))
    } else {
      data_input_cd = conv_to_cd_data(data_input) # convert data for CD algo
      fit.cd = cd.run(indata = data_input_cd)
    }

    A_est_by_lam_cd = calc_A_est_by_lam_cd(fit.cd)
    temp_cd_result = eval_by_lam(A_est_by_lam_cd, A_true, n_edges)
    temp_cd_sel = mod_sel(temp_cd_result, n_edges, sel_crit = sel_crit)
    result_cd = rbind(result_cd, temp_cd_sel)

  }

  total_result_n_edges[[graph_type]] <- n_edges
  total_result_cd[[graph_type]] <-  result_cd
  total_result_noteargis[[graph_type]] <- result_noteargis

}

sink("table3.txt")
for (g in 1:n_types){ # g = 1
  graph_g = graph_types[g]
  cat("graph:", graph_g)
  cat("\n")
  cd_ave <- as.matrix(t(round(colMeans(total_result_cd[[graph_g]]),2)))
  #rownames(cd_ave) = c("CD")

  noteargis_ave <-  as.matrix(t(round(
    colMeans(total_result_noteargis[[graph_g]]),2)))
  #rownames(noteargis_ave) = "NOTEARGIS"

  final_result = rbind(noteargis_ave, cd_ave)
  rownames(final_result) <- c("NOTEARGIS", "    CD")
  print(final_result)
  cat("========================================================\n")
}
sink()

# _ (3) Generate the Table 4 ####

use_saved_data = 1
# if 1, it will used the saved results to generate table.
# if 0, it will run the functions (it will take much longer)

n_obs = 200
n_nodes = 30
types_by_node = rep("m", n_nodes)

iterations = 50
seed_para = 1
graph_types = c("bi","rand", "sf", "sw","tree")
methods = c("noteargis", "cd")
n_types = length(graph_types)
n_methods = length(methods)

n_resps_by_node = gen_node_levels(types_by_node, max_n_levels = 4, seed = seed_para)

total_result_cd = as.list(rep(0,n_types))
names(total_result_cd) <- graph_types
total_result_noteargis <- total_result_cd
total_result_n_edges <- total_result_cd

sel_crit = "JI"

for (graph_type in graph_types){ # graph_type = "bi"

  A_true = gen_adj_mat(n_nodes, graph_type, seed = seed_para)
  n_edges = sum(A_true)

  W_true = gen_para(A_true, types_by_node, n_resps_by_node, seed = seed_para)
  result_cd = data.frame()
  result_noteargis = data.frame()

  for (seed_data in 1:iterations){ # seed_data = 1

    ## NOTEARGIS
    if (use_saved_data == 1){
      obj_name = paste0("noteargis_", graph_type, "_data_",seed_data,"_table4.rds")
      fit.noteargis = readRDS(system.file("extdata", obj_name,
                                          package = "noteargisEval"))
    } else {
      data_input = gen_data(n_obs, A_true, W_true, seed = seed_data)
      fit.noteargis = noteargis(data_input, "m", n_resps_by_node, verbose = T)
    }

    temp_noteargis_result = eval_by_lam(fit.noteargis$A_est_by_lam, A_true, n_edges)
    temp_noteargis_sel = mod_sel(temp_noteargis_result, n_edges, sel_crit = sel_crit)
    result_noteargis = rbind(result_noteargis, temp_noteargis_sel)

    ## CD
    if (use_saved_data == 1){
      obj_name = paste0("cd_", graph_type, "_data_",seed_data,"_table4.rds")
      fit.cd = readRDS(system.file("extdata", obj_name, package = "noteargisEval"))
    } else {
      data_input_cd = conv_to_cd_data(data_input) # convert data for CD algo
      fit.cd = cd.run(indata = data_input_cd)
    }

    A_est_by_lam_cd = calc_A_est_by_lam_cd(fit.cd)
    temp_cd_result = eval_by_lam(A_est_by_lam_cd, A_true, n_edges)
    temp_cd_sel = mod_sel(temp_cd_result, n_edges, sel_crit = sel_crit)
    result_cd = rbind(result_cd, temp_cd_sel)

  }

  total_result_n_edges[[graph_type]] <- n_edges
  total_result_cd[[graph_type]] <-  result_cd
  total_result_noteargis[[graph_type]] <- result_noteargis

}

sink("table4.txt")
for (g in 1:n_types){ # g = 1
  graph_g = graph_types[g]
  cat("graph:", graph_g)
  cat("\n")
  cd_ave <- as.matrix(t(round(colMeans(total_result_cd[[graph_g]]),2)))
  #rownames(cd_ave) = c("CD")

  noteargis_ave <-  as.matrix(t(round(
    colMeans(total_result_noteargis[[graph_g]]),2)))
  #rownames(noteargis_ave) = "NOTEARGIS"

  final_result = rbind(noteargis_ave, cd_ave)
  rownames(final_result) <- c("NOTEARGIS", "    CD")
  print(final_result)
  cat("========================================================\n")
}
sink()

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
## Section 5.5 Application: Alarm Data ####
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

use_saved_data = 1
# if 1, it will used the saved results to generate table.
# if 0, it will run the functions (it will take much longer)

# _ (1) Generate the Table 5 ####
data(alarm)
data(alarm_true_str)


n_obs = 200 # number of susbset
n_nodes = ncol(alarm)
types_by_node = rep("m", n_nodes)

data_raw = alarm
n_raw = nrow(data_raw)
n_levels_raw = sapply(data_raw, function(x) nlevels(x))

iterations = 20
data_seeds_init = 1:iterations
data_seeds = c()

## find out subsetting seeds to guanrantee that a subset also have the same
##  ... number of levels as n_levels_raw
for (i in data_seeds_init){

  seed_data <- i
  kk = 0
  while (TRUE){
    set.seed(seed_data)
    sel = sample(1:n_raw, n_obs)
    data_input = data_raw[sel,]
    n_levels = sapply(data_input, function(x) length(unique(x)))
    if(all(n_levels_raw == n_levels)) break
    kk = kk + 1
    seed_data = iterations * kk + seed_data
  }
  data_seeds[i] <- seed_data
}

graph_types = c("alarm")
methods = c("noteargis", "cd")
n_types = length(graph_types)
n_methods = length(methods)

n_resps_by_node = n_levels_raw

total_result_cd = as.list(rep(0,n_types))
names(total_result_cd) <- graph_types
total_result_noteargis <- total_result_cd
total_result_n_edges <- total_result_cd

sel_crit = "JI"

for (graph_type in graph_types){ # graph_type = "alarm"

  A_true = alarm_true_str$A_true
  n_edges = alarm_true_str$n_edges

  result_cd = data.frame()
  result_noteargis = data.frame()

  for (seed_data in data_seeds){ # seed_data = 1

    ## NOTEARGIS
    if (use_saved_data == 1){
      obj_name = paste0("noteargis_", graph_type, "_data_",seed_data,"_table5.rds")
      fit.noteargis = readRDS(system.file("extdata", obj_name,
                                          package = "noteargisEval"))
    } else {
      set.seed(seed_data)
      sel = sample(1:n_raw, n_obs)
      data_input = data_raw[sel,]
      fit.noteargis = noteargis(data_input, "m", n_resps_by_node, verbose = T)
    }

    temp_noteargis_result = eval_by_lam(fit.noteargis$A_est_by_lam, A_true, n_edges)
    temp_noteargis_sel = mod_sel(temp_noteargis_result, n_edges, sel_crit = sel_crit)
    result_noteargis = rbind(result_noteargis, temp_noteargis_sel)

    ## CD
    if (use_saved_data == 1){
      obj_name = paste0("cd_", graph_type, "_data_",seed_data,"_table5.rds")
      fit.cd = readRDS(system.file("extdata", obj_name, package = "noteargisEval"))
    } else {
      data_input_cd = conv_to_cd_data(data_input) # convert data for CD algo
      fit.cd = cd.run(indata = data_input_cd)
    }

    A_est_by_lam_cd = calc_A_est_by_lam_cd(fit.cd)
    temp_cd_result = eval_by_lam(A_est_by_lam_cd, A_true, n_edges)
    temp_cd_sel = mod_sel(temp_cd_result, n_edges, sel_crit = sel_crit)
    result_cd = rbind(result_cd, temp_cd_sel)

  }

  total_result_n_edges[[graph_type]] <- n_edges
  total_result_cd[[graph_type]] <-  result_cd
  total_result_noteargis[[graph_type]] <- result_noteargis

}

sink("table5.txt")
for (g in 1:n_types){ # g = 1
  graph_g = graph_types[g]
  cat("graph:", graph_g)
  cat("\n")
  cd_ave <- as.matrix(t(round(colMeans(total_result_cd[[graph_g]]),2)))

  noteargis_ave <-  as.matrix(t(round(
    colMeans(total_result_noteargis[[graph_g]]),2)))

  final_result = rbind(noteargis_ave, cd_ave)
  rownames(final_result) <- c("NOTEARGIS", "    CD")
  print(final_result)
  cat("========================================================\n")
}
sink()

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
## Appendix D.1 Generate Graph Structure ####
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
n_nodes = 4
(A_true = gen_adj_mat(n_nodes, "rand", seed = 1))

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
## Appendix D.2 Generate Parameters ####
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
types_by_node = rep("c", n_nodes) # grouped continuous variables
n_resps_by_node = c(2,1,3,2) # number of elements for each group
(W_true_grp_conti = gen_para(A_true, types_by_node, n_resps_by_node))

types_by_node = rep("m", n_nodes) # multinomial responses
n_resps_by_node = c(2,4,2,3) # number of levels
W_true_multi = gen_para(A_true, types_by_node, n_resps_by_node)
round(W_true_multi,2)

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
## Appendix D.3 Generate Dataset ####
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
(data_grp_conti = gen_data(n_obs = 5, A_true, W_true_grp_conti, seed = 1))
(data_multi = gen_data(n_obs = 5, A_true, W_true_multi, seed = 2))
