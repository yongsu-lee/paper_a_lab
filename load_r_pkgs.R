## Install bioconductor package first, 
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install()

# BiocManager::install(c("graph", "RBGL"))


# if(!require(pcalg)){
#   install.packages("pcalg", quiet = T)
#   library("pcalg", quietly = T)
# }

##

if(!require(VGAM)){
  install.packages("VGAM", quiet = T)
  library("VGAM", quietly = T)
}


if(!require(Matrix)){
  install.packages("Matrix", quiet = T)
  library("Matrix", quietly = T)
}

if(!require(igraph)){
  install.packages("igraph", quiet = T)
  library("igraph", quietly = T)
}

if(!require(Rcpp)){
  install.packages("Rcpp", quiet = T)
  library("Rcpp", quietly = T)
}

if(!require(RcppArmadillo)){
  install.packages("RcppArmadillo", quiet = T)
  library("RcppArmadillo", quietly = T)
}

# if(!require(RcppEigen)){
#   install.packages("RcppEigen", quiet = T)
#   library("RcppEigen", quietly = T)
# }


if(!require(lbfgs)){
  install.packages("lbfgs", quiet = T)
  library("lbfgs", quietly = T)
}

# if(!require(nloptr)){
#   install.packages("nloptr", quiet = T)
#   library("nloptr", quietly = T)
# }

if(!require(expm)){
  install.packages("expm", quiet = T)
  library("expm", quietly = T)
}


## For comparison purposes
if(!require(deal)){
  install.packages("deal", quiet = T)
  library("deal", quietly = T)
}

if(!require(nnet)){
  install.packages("nnet", quiet = T)
  library("nnet", quietly = T)
}


## For load_cd_realted_cdoes.R ----

if(!require(discretecdAlgorithm)){
  install.packages("discretecdAlgorithm", quiet = T)
  library("discretecdAlgorithm", quietly = T)
}

if(!require(sparsebn)){
  install.packages("sparsebn", quiet = T)
  library("sparsebn", quietly = T)
}

## Imputing missing values

if(!require(bnstruct)){
  install.packages("bnstruct", quiet = T)
  library("bnstruct", quietly = T)
}

