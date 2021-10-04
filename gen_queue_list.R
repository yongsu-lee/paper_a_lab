queue_list = function(args_list){
  
  n_args = length(args_list)
  combn = expand.grid(args_list)
  write.table(combn, file = "queue_list_simu1_more_iter", sep=", ", 
              quote = F, row.names = F, col.names = F)
}

args_list = list(c("bi", "rand", "sf", "sw", "tree"), 51:100)


# args_list = list(1:20, c("alarm"), c("FALSE"))
queue_list(args_list)