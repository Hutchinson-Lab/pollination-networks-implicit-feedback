# ============================================================================ #
# Evaluation metrics
#  1. MPR
#  2. Precision-recall (PR)
# Author: Eugene Seo 
# Data: February 24, 2018
# ============================================================================ #

EvaluateCommon <- function(Ctrain, test.mat, new.mat, Atest, pred.mat, target) {
  dimnames(pred.mat) <- dimnames(Ctrain)  
  pred.mat <- GetCommonMatrics(pred.mat, test.mat)[[1]] 
  int.mat <- pred.mat * Atest 
  
  if(target == "All") {
    res <- EvaluationAll(test.mat, int.mat, new.mat, "All")
  } else if(target == "New") {
    res <- EvaluationAll(test.mat, int.mat, new.mat, "New")
  } 
  return(res)
}

ValidateCommon <- function(Ctrain, Ctest, Cnew, Atest, pred.mat, target) {
  dimnames(pred.mat) <- dimnames(Ctrain)  
  pred.mat <- GetCommonMatrics(pred.mat, Ctest)[[1]] 
  int.mat <- pred.mat * Atest 
  
  if(target=="All") {  
    res <- ComputeMPR(Ctest, int.mat, Cnew, "All")
  } else if(target=="New") {  
    res <- ComputeMPR(Ctest, int.mat, Cnew, "New")
  } 
  return(res)
}

EvaluationAll <- function(test.mat, int.mat, new.mat, target) {
  intNum <- rowSums((test.mat>0)*1)
  remove_idx <- which(intNum==0)
  if(length(remove_idx) > 0) {
    test.mat <- test.mat[-remove_idx, ]
    int.mat <- int.mat[-remove_idx, ]
    new.mat <- new.mat[-remove_idx, ]
  }  
  mpr <- ComputeMPR(test.mat, int.mat, new.mat, target)
  pr <- ComputePR(test.mat, int.mat, new.mat, target)
  return(list(MPR=mpr, PR=pr))
}

ComputeMPR <- function(test.mat, int.mat, new.mat, target) {
  if(target == "All") {    
    rank.mat <- RankByRow(int.mat, "average")
    perc.rank <- PercentileRanking(rank.mat)
  } else if(target == "New") {    
    pos.idx = which(test.mat>0)
    new.idx = which(new.mat>0)
    train.idx = pos.idx[!pos.idx %in% new.idx]
    int.mat[train.idx] = 0
    
    intNum <- rowSums((new.mat>0)*1) # pollinators having no new interactions
    remove_idx <- which(intNum==0)
    if(length(remove_idx) > 0) {
      test.mat <- test.mat[-remove_idx, ]
      int.mat <- int.mat[-remove_idx, ]
      new.mat <- new.mat[-remove_idx, ]
    }
    rank.mat <- RankByRow(int.mat, "average")
    perc.rank <- PercentileRanking(rank.mat)
  }
  mpr <- MPR(test.mat, perc.rank, new.mat, target)
  return(mpr)
}

ComputePR <- function(test.mat, pred.mat, new.mat, mode) {  
  topN.max <- ncol(test.mat)
  PR.data <- data.frame(matrix(NA, ncol=2))
  
  if(mode == "All") {    
    rank.mat <- RankByRow(pred.mat, "random")
    relevant.num <- length(which(test.mat>0))
  } else if(mode == "New") {   
    pos.idx <- which(test.mat>0)
    new.idx <- which(new.mat>0)
    train.idx <- pos.idx[!pos.idx %in% new.idx]
    pred.mat[train.idx] <- 0
    
    intNum <- rowSums((new.mat>0)*1)
    remove_idx <- which(intNum==0)
    if(length(remove_idx) > 0) {
      test.mat <- test.mat[-remove_idx, ]
      pred.mat <- pred.mat[-remove_idx, ]
      new.mat <- new.mat[-remove_idx, ]
    }
  
    rank.mat <- RankByRow(pred.mat, "random")
    relevant.num <- length(which(new.mat>0))
  }
  
  for(topN in 1:topN.max) {
    recommend.num <- 0
    correct.num <- 0      
    idx <- which(rank.mat <= topN)
    recommend.num <- recommend.num + length(idx)
    if(mode == "All") {
      correct.num <- correct.num + length(which(test.mat[idx]>0))
    } else if(mode == "New") {
      correct.num <- correct.num + length(which(new.mat[idx]>0))
    }  
    Precision <- correct.num / recommend.num
    Recall <- correct.num / relevant.num
    PR.data[topN, 1] <- Precision
    PR.data[topN, 2] <- Recall
  }  
  F1 <- (1.25 * PR.data[,1] * PR.data[,2]) / (0.25 * PR.data[,1] + PR.data[,2])
  PR.data <- cbind(PR.data, F1)
  colnames(PR.data) <- c("Precision", "Recall", "F1") 
  return(PR.data)
}

RankByRow <- function(mat, method) {
  rank.mat <- t(apply(-mat, 1, rank, ties.method=method))
  return(rank.mat)
}

PercentileRanking <- function(x) {
  amin <- min(x)
  amax <- max(x)
  x <- (x-amin) / (amax - amin) * 100
  return(x)
}

MPR <- function(test.mat, perc.rank, new.mat, target) {
  if(target == "All") {
    mpr <- sum(test.mat * perc.rank)/sum(test.mat)
  } else if(target =="New") {
    idx <- which(new.mat>0)
    mpr <- sum(test.mat[idx] * perc.rank[idx]) / sum(test.mat[idx])
  }
  return(mpr)
}