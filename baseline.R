# ============================================================================ #
# Baseline methods
#  1. availability: The prefrence is determined by flower availability 
#  2. popularity: The prefrence is determined by flower popularity       
#  3. ItemkNN: The prefrence is determined by the similarity of items
#  4. UserkNN: The prefrence is determined by the similatiy of users                                                     
# Author: Eugene Seo 
# Data: February 24, 2018
# ============================================================================ #
source("knn.R")

Baseline <- function(train.mat, op, flw.test.abund=0) {
  if (op == "availability") {
    pred.mat <- t(matrix(rep(flw.test.abund, nrow(train.mat)), 
                             ncol=nrow(train.mat),))  
  } else if (op == "popularity") {
    pred.mat <- t(replicate(nrow(train.mat), colSums(train.mat)))
  } else if (op == "ItemNN") {
    neighbors <- findKNN(train.mat, "item") 
    pred.mat <- prediction(train.mat, neighbors$kNN, neighbors$sim, 
                           ncol(train.mat)-1, TRUE, "item")
  } else if (op == "UserNN") {
    neighbors <- findKNN(train.mat, "user") 
    pred.mat <- prediction(train.mat, neighbors$kNN, neighbors$sim, 
                           nrow(train.mat)-1, TRUE, "user")
  }
  return(pred.mat)
}