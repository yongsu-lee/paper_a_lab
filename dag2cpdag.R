## Function from the package pcalg
dag2cpdag <- function(adj)
{
  ## Purpose: Compute the (unique) completed partially directed graph (CPDAG)
  ## that corresponds to the input DAG; result is a graph object
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## - dag: input DAG (graph object)
  ## ----------------------------------------------------------------------
  ## Author: Diego Colombo, Date: 10 Jun 2013, 11:06
  # adj = A_true
  
  amat <- adj
  amat[amat != 0] <- 1
  skel.amat <- amat + t(amat)
  skel.amat[skel.amat == 2] <- 1
  cpdag <- skel.amat
  
  ## search the v-structures in the DAG
  ind <- which((amat == 1 & t(amat) == 0), arr.ind = TRUE)
  tripleMatrix <- matrix(,0,3)
  ## Go through all edges
  for (i in seq_len(nrow(ind))) { ## MM(FIXME): growth of tripleMatrix
    x <- ind[i,1]
    y <- ind[i,2]
    indY <- setdiff(which((amat[,y] == 1 & amat[y,] == 0), arr.ind = TRUE),x) ## x-> y <- z
    if(length(newZ <- indY[amat[x,indY] == 0])) ## deparse.l.=0: no colnames
      tripleMatrix <- rbind(tripleMatrix, cbind(x, y, newZ, deparse.level=0),
                            deparse.level=0)
  }
  if ((m <- nrow(tripleMatrix)) > 0) {
    deleteDupl <- logical(m)# all FALSE
    for (i in seq_len(m))
      if (tripleMatrix[i,1] > tripleMatrix[i,3])
        deleteDupl[i] <- TRUE
      if(any(deleteDupl))
        tripleMatrix <- tripleMatrix[!deleteDupl,, drop=FALSE]
      
      ## orient the v-structures in the CPDAG
      for (i in seq_len(nrow(tripleMatrix))) {
        x <- tripleMatrix[i,1]
        y <- tripleMatrix[i,2]
        z <- tripleMatrix[i,3]
        cpdag[x,y] <- cpdag[z,y] <- 1
        cpdag[y,x] <- cpdag[y,z] <- 0
      }
  }
  
  ## orient the edges with the 3 orientation rules
  repeat {
    old_cpdag <- cpdag
    ## Rule 1
    ind <- which((cpdag == 1 & t(cpdag) == 0), arr.ind = TRUE)
    for (i in seq_len(nrow(ind))) {
      a <- ind[i, 1]
      b <- ind[i, 2]
      isC <- ((cpdag[b, ] == 1 & cpdag[, b] == 1) &
                (cpdag[a, ] == 0 & cpdag[, a] == 0))
      if (any(isC)) {
        indC <- which(isC)
        cpdag[b, indC] <- 1
        cpdag[indC, b] <- 0
      }
    }
    ## Rule 2
    ind <- which((cpdag == 1 & t(cpdag) == 1), arr.ind = TRUE)
    for (i in seq_len(nrow(ind))) {
      a <- ind[i, 1]
      b <- ind[i, 2]
      isC <- ((cpdag[a, ] == 1 & cpdag[, a] == 0) &
                (cpdag[, b] == 1 & cpdag[b, ] == 0))
      if (any(isC)) {
        cpdag[a, b] <- 1
        cpdag[b, a] <- 0
      }
    }
    ## Rule 3
    ind <- which((cpdag == 1 & t(cpdag) == 1), arr.ind = TRUE)
    for (i in seq_len(nrow(ind))) {
      a <- ind[i, 1]
      b <- ind[i, 2]
      indC <- which((cpdag[a, ] == 1 & cpdag[, a] == 1) &
                      (cpdag[, b] == 1 & cpdag[b, ] == 0))
      if (length(indC) >= 2) {
        cmb.C <- combn(indC, 2)
        cC1 <- cmb.C[1, ]
        cC2 <- cmb.C[2, ]
        for (j in seq_along(cC1)) {
          c1 <- cC1[j]
          c2 <- cC2[j]
          if (c1 != c2 && cpdag[c1, c2] == 0 && cpdag[c2,c1] == 0) {
            cpdag[a, b] <- 1
            cpdag[b, a] <- 0
            break
          }
        }
      }
    }
    if (all(cpdag == old_cpdag))
      break
  }
  
  cpdag_adj = cpdag
  
  return(cpdag_adj)
  
}