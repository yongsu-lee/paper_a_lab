rm(list = ls()); 
source("pkg_ctrl.R")
source("sourcing_nocpp.R")
sourceCpp("./lkhd_mix.cpp")

data_input = readRDS("data_input.rds")

n_lams = 30

graph_est_by_lam = as.list(rep(NA,n_lams))
pushed_A_est_by_lam = as.list(rep(NA,n_lams))
pushed_A_est_by_lam = as.list(rep(NA,n_lams))
attr(pushed_A_est_by_lam, "n_edges") <- c()
pushed_Beta_est_by_lam = pushed_A_est_by_lam

temp = readRDS(paste0("./app_result2/result_path_", 1, ".rds"))
data_info = temp

## load server result
for (ell in 1:n_lams){ # ell = 1
  
  file_name = paste0("./app_result2/result_path_", ell, ".rds")
  result = readRDS(file_name)
  
  graph_est_by_lam[[ell]] <- result$graph_est_by_lam[[1]]
  pushed_A_est_by_lam[[ell]] <- result$pushed_A_est_by_lam[[1]]
  attr(pushed_A_est_by_lam, "n_edges")[ell] <-  
    attr(result$pushed_A_est_by_lam, "n_edges")
  
  pushed_Beta_est_by_lam[[ell]] <- result$pushed_Beta_est_by_lam[[1]]
  
}

tun_sel_lkhd(pushed_A_est_by_lam, pushed_Beta_est_by_lam, data_info)

var_names_temp = names(data_input)

## delete alcohol and smoking
var_names = var_names_temp[-c(3,10)] 

var_names_short <- c("gend", "symp", "hbsag", "hbeag", "hbcab", "hcvab", "cirr",
                     "ende", "diab", "obes", "hech", "aht", "cri", "hiv", "nash", "vari",
                     "spln", "pht", "pvt", "meta", "hmrk", "age", "alch", "cig",
                     "ps", "enc", "asc", "inr", "afp", "hemo", "mcv", "leuc",
                     "plat", "albu", "tot.bil", "alt", "ast", "ggt", "alp", "tp",
                     "crea", "nodules", "major", "dir.bil", "iron", "sat", "ferr", "class")

cbind(var_names, var_names_short)

## original graph 
pdf("app_original.pdf")
par(mar=c(0,0,0,0))
plot(graph_est_by_lam[[17]], 
     vertex.label = var_names_short, 
     vertex.label.cex=rep(1,50),
     vertex.size = 10,
     edge.arrow.size=0.3,
     )
dev.off()



## excluding isolated nodes (no kid and no parent)
A_est_reduced = pushed_A_est_by_lam[[17]]
combined = cbind(colSums(A_est_reduced), rowSums(A_est_reduced))
sig_nodes = rowSums(combined) !=0

# sig_nodes = c(11,30,29,15,24,10,13,21,18,20,2,3,9,22,25,16)

A_est_reduced = A_est_reduced[sig_nodes, sig_nodes]
graph_est_reduced = graph_from_adjacency_matrix(A_est_reduced)

var_names_reduced = var_names[sig_nodes]

plot(graph_est_reduced , 
     vertex.label = var_names_reduced , 
     vertex.label.cex=rep(0.8,50),
     vertex.size = 10,
     edge.arrow.size=0.3
)


which(A_est_reduced[,29] != 0)

idx_roots = which(colSums(A_est_reduced) == 0)
var_names_reduced[idx_roots]
idx_terminals = which(rowSums(A_est_reduced) == 0)
var_names_reduced[idx_terminals]

var_names_reduced[25]



sig_nodes = c(11,30,29,15,24,10,13,21,18,20,2,3,9,22,25,16)

A_est_reduced = A_est_reduced[sig_nodes, sig_nodes]
graph_est_reduced = graph_from_adjacency_matrix(A_est_reduced)

var_names_reduced = var_names_reduced[sig_nodes]

plot(graph_est_reduced , 
     vertex.label = var_names_reduced , 
     vertex.label.cex=rep(0.8,50),
     vertex.size = 10,
     edge.arrow.size=0.5
)


which(A_est_reduced[,2] != 0)

idx_roots = which(colSums(A_est_reduced) == 0)
var_names_reduced[idx_roots]
idx_terminals = which(rowSums(A_est_reduced) == 0)
var_names_reduced[idx_terminals]

var_names_reduced[25]

