# ============================================================================ #
# kNN method
#  1. Item-based: Find the similar items to predict the prefernce score. 
#  2. User-based: Find the similar users to predict the prefernce score.                                     
# Author: Eugene Seo 
# Data: February 24, 2018
# ============================================================================ #

library("proxy")
findKNN <- function(mat, mode) { 
  if(mode=="user") {
    mat <- t(mat)
  }
  sim.mat <- as.matrix(simil(mat, method="Cosine", diag=TRUE, upper=TRUE, 
                             by_rows=FALSE), diag=1)
  n <- nrow(sim.mat)
  k <- ncol(sim.mat)-1
  knn.mat <- matrix(0, nrow=n, ncol=k)
  for(i in 1:n) {
    neighbors <- order(sim.mat[i, ], decreasing=TRUE)
    neighbors <- neighbors[neighbors != i]
    knn.mat[i,] <- neighbors[1:k]    
  }
  return(list(kNN=knn.mat, sim=sim.mat))
}

prediction <- function(mat, knn.mat, sim.mat, k, onlyrated, mode) {  
  pred.mat <- matrix(0, nrow=nrow(mat), ncol=ncol(mat))  
  if(mode == "item") {
    if(k >= ncol(mat)) {
      print("k is large than the number of items.")
      return(FALSE)
    }
    for(i in 1:nrow(mat)) {
      rated.items <- which(mat[i, ]>0)
      for(j in 1:ncol(mat)) {      
        similar.items <- knn.mat[j, 1:k]
        if(onlyrated) {
          common.items <- intersect(rated.items, similar.items)
        } else {
          common.items <- similar.items
        }
        similarity <- sim.mat[j, common.items]
        user.ratings <- mat[i, common.items]
        if(sum(similarity) == 0) {
          pred.mat[i,j] <- 0
        } else {
          pred.mat[i,j] <- sum(similarity * user.ratings) / sum(similarity)
        }
      }
    }
  } else if (mode == "user") {
    if(k >= nrow(mat)) {
      print("k is large than the number of user.")      
      return(FALSE)
    }
    for(i in 1:nrow(mat)) {
      similar.users <- knn.mat[i, 1:k]
      for(j in 1:ncol(mat)) {
        rated.users <- which(mat[, j]>0)
        if(onlyrated) {
          common.users <- intersect(rated.users, similar.users)
        } else {
          common.users <- similar.users
        }
        similarity <- sim.mat[i, common.users]
        common.users.ratings <- mat[common.users, j]
        if(sum(similarity) == 0) {
          pred.mat[i,j] <- 0
        } else {
          pred.mat[i,j] <- sum(similarity * common.users.ratings) / 
                           sum(similarity)
        }
      }
    }
  }
  rownames(pred.mat) <- rownames(mat)
  colnames(pred.mat) <- colnames(mat)
  return(pred.mat)
}