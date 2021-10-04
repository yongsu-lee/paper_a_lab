## Load necessary packages
source("./load_r_pkgs.R")
#.. More details of the usage of each package can be found in "pkg_ctrl.R" file.

## Load commonly-used functions (short-length functions)
source("./shared_fun.R")

## Cpp-enhanced functions to deal with bottlenecks
sourceCpp("./arma_trexpm.cpp")
sourceCpp("./obj0_mix.cpp")
sourceCpp("./grad0_mix.cpp")
sourceCpp("./lkhd_mix.cpp")
sourceCpp("./lkhd_multi.cpp")
sourceCpp("./lkhd_ordin.cpp")
sourceCpp("./lkhd_conti.cpp")

## Generating data source
source("./gen_graph_adj_mat.R")
source("./gen_para.R")
source("./gen_data.R")

## Main algorithm
source("./glmdag.R")
source("./initialize_W.R")
source("./gen_lambdas.R")
source("./admm_loop.R")
source("./mix.W_update.R")
source("./Beta_update.R")

## Evaluation source
source("./dag2cpdag.R")
source("./push_dag.R")
source("./tun_sel_lkhd.R")
# source("./eval_metrics.R")


## Load old function for comparison








