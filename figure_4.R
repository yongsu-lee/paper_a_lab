#######################################################
## Figure 4 (Large Graphs, Nom vs Ord) ####
#######################################################

seed_para = 1

n_conti = 16
n_multi = 14
n_nodes = n_conti + n_multi 

graph_types = c("bi", "rand", "sf", "sw","tree")
graph_names = c("Bipartite", "Random DAG", "Scale-Free", "Small-World", "Tree")

types_by_node = gen_node_types(n_conti, n_multi, 0, seed = seed_para)

m = c(1,1,2,2,3,3,1,1,2,2,3,3,0,4,4,5,5,0,0,4,4,5,5,0)
M = matrix(m, nrow =4, ncol = 6, byrow = T); layout(mat = M)

pdf("./results/Figure_4.pdf", width = 10)
layout(mat = M); par(mar=c(0,0.5,3,0.5)); par(oma=c(0,0,0,0))

for (j in 1:5){
  
  graph_type = graph_types[j]
  graph_name = graph_names[j]
  
  graph_set = gen_graph_adj_mat(n_nodes, graph_type, seed = seed_para)
  graph_true = graph_set$graph_true
  A_true = graph_set$A_true
  
  fancy_dag_plot(A_true, types_by_node, title = graph_name, legend = F, 
                 edge.arrow.size = 0.3,
                 vertex.size = 10,
                 layout = layout_set(graph_true, graph_type))
  
}

dev.off()
