#######################################################
## Figure 5 (Original Est Graph with HCC) ####
#######################################################

if(T){ # for debugging in local, delete when submitting
  
  setwd("~/iCloud/paper_a/paper_a_submit/biom_j/biom_j_code_submit/")
  rm(list=ls())
  source("./codes/functions/00_load_ftns.R")
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  
  queue_args = read.table("./codes/functions/queue_list_hcc", sep =",", 
                         strip.white = T)
  ell = queue_args[1,]
  n_lams = 30 # the number of tuning parameters
  save_file_name = paste0("./results/hcc/hcc_path_",ell,".rds")
  
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
  ell = as.integer(queue_args[1])
  save_file_name = paste0("hcc_path_",ell,".rds")

  n_lams = 30
  
}

## Import hcc data

# Get hcc dataset from UCI ML Repository
url = "https://archive.ics.uci.edu/ml/machine-learning-databases/00423/hcc-survival.zip"
download.file(url, destfile = file.path(".", basename(url)))
unzip(zipfile = "hcc-survival.zip", exdir = ".")
file.copy("./hcc-survival/hcc-data.txt", ".")

# read original data
hcc_temp <- read.table("hcc-data.txt", sep=",", na.strings = "?")
n = nrow(hcc_temp)
hcc_temp_b4_impute = as.matrix(hcc_temp)

# Export the original data with xlsx format
hcc_info <- read.csv("hcc_info.csv") 
names(hcc_temp) <- hcc_info$var_names


## Data cleaning
source("hcc_cleaning.R")
#.. See details in "./codes/hcc_cleaning.R"
data_input = hcc_imputed
# write.csv(hcc_imputed, "hcc_imputed.csv")

## Run GLMDAG
terminal_nodes = c("class")
root_nodes = c("gend", "ende")

result = glmdag(data_input, n_lams = n_lams, 
                root_nodes = root_nodes, terminal_nodes = terminal_nodes,
                path_par = T, path_par_num = ell, verbose = T)

saveRDS(result, file = save_file_name)




