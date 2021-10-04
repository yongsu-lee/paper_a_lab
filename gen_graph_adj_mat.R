#' Generate a true graph and an adjacency matrix
#'
#' @description Generate a desired graph (with \code{igraph} format) with corresponding adjacency matrix as an binary matrix. If \eqn{(i,j)} is one, it implies there exists an edge \eqn{i -> j}.
#'
#' @param n_nodes positive integer scalar (larger than 1), the number of nodes of the desired graph.
#' @param graph_type character scalar, one of types of graphs in the simulation studies. See details below for more information.
#' @param seed integer scalar, in order to reproduce the same graph.
#'
#' @return \item{graph_true}{a requested \code{igraph} graph} \item{A_true}{a corresponding adjacecny matrix of the \code{graph_true}}
#'
#' @export
#' @import igraph pcalg
#' @details There are seven graph types in the simulation studies: "bi" for Bipartite graph, "sf" for Scale-Free graph, "star" for Star-Out graph, "sw" for Small-World, "tree" for the tree with two childeren, "rand" for random DAG, and "chain" for a Chain graph. Except "rand" and "chain", all the other graphs generated based on R package \strong{igraph}. The Random DAG is generaated based on R package \strong{pcalg}.
#' @examples
#' gen_graph_adj_mat(n_nodes = 10, graph_type = "bi")
gen_graph_adj_mat = function(
  n_nodes, graph_type=c("bi","sf","star","sw","tree","rand","chain"), seed=NA)
{

  graph_type = match.arg(graph_type)
  if (is.na(seed)) seed = 1
  
  if (n_nodes == 2){
    A_true = matrix(c(0,0,1,0), 2, 2)
    graph_true <- igraph::graph_from_adjacency_matrix( A_true )
    A_true = igraph::as_adjacency_matrix(graph_true)
    
  } else {
    
    while(TRUE){
      
      switch(graph_type,
             
             "bi" = {
               if (!is.na(seed)) set.seed(seed)
               graph_true = igraph::sample_bipartite(
                 n1 = round(0.5*n_nodes), n2 = round(0.5*n_nodes), 
                 type = "Gnm", m = n_nodes, directed = TRUE, mode = "in")
               A_true = igraph::as_adjacency_matrix(graph_true) },
             
             "tree" = {
               if (!is.na(seed)) set.seed(seed)
               graph_true = igraph::graph.tree(n = n_nodes, children = 2)
               A_true = igraph::as_adjacency_matrix(graph_true) },
             
             "sf" = {
               if (!is.na(seed)) set.seed(seed)
               graph_true = igraph::sample_pa(n = n_nodes, directed = T)
               A_true = igraph::as_adjacency_matrix(graph_true) },
             
             "sw" = {
               if (!is.na(seed)) set.seed(seed)
               graph_temp <- igraph::sample_smallworld(dim = 1, size = n_nodes, 
                                                       nei = 2,  p = 0.5)
               A_temp <- igraph::as_adjacency_matrix(graph_temp)
               graph_true = igraph::as.directed(graph_temp, mode="arbitrary")
               A_true = igraph::as_adjacency_matrix(graph_true)},
             
             "rand" = {
               if (!is.na(seed)) set.seed(seed)
               graph_true = erdos.renyi.game(n_nodes, n_nodes, 
                                             type = "gnm", directed = T)
               A_true = igraph::as_adjacency_matrix(graph_true) },
             
             "star" = {
               if (!is.na(seed)) set.seed(seed)
               graph_true = igraph::make_star(n = n_nodes, mode = "out")
               A_true = igraph::as_adjacency_matrix(graph_true)
               A_true = as.matrix(A_true) },
             
             "chain" = {
               if (!is.na(seed)) set.seed(seed)
               A_true = matrix(0, n_nodes, n_nodes)
               A_true[(row(A_true) + 1) == col(A_true) ] <-1
               graph_true = igraph::graph_from_adjacency_matrix(A_true) }
      )
      
      if (igraph::is_dag(graph_true)) break
      
      seed <- seed + 1
      
    }
  }
  
  attr(graph_true, "graph_type") = graph_type
  A_true = as.matrix(A_true)
  colnames(A_true) = paste0("Z",1:n_nodes)
  rownames(A_true) = paste0("Z",1:n_nodes)
  
  func_result = list(graph_true = graph_true, A_true = A_true, seed = seed)
  
  return(func_result)
  
}
