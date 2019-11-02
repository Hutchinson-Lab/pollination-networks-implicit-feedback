# ============================================================================ #
# Separate data into a validation set and a test set
# (e.g. test set - 2015
#       validation set - 2014
#       training set for validation - 2011, 2012, 2013       
#       training set for testing - 2011, 2012, 2013, 2014
# Author: Eugene Seo 
# Data: February 24, 2018
# ============================================================================ #
source("build_flower_avail.R")

GetValidEval <- function(test.y) {
  # Set the validation set
  if(test.y == year.list[1]){
    valid.test.y <- tail(year.list, n=1)
  } else {
    valid.test.y <- test.y - 1
  }  
  # Set the training set for validation
  valid.train.y <- year.list[!year.list== valid.test.y & !year.list== test.y]
  # Set the training set for testing
  train.y <- year.list[!year.list==test.y]
  
  cat(c("For valiation,\n"))
  cat(c("\ttrain.years:", valid.train.y, "\n\ttest.year:", valid.test.y, "\n"))
  cat(c("For evaluation,\n"))
  cat(c("\ttrain.years:", train.y, "\n\ttest.year:", test.y, "\n\n"))
  
  cat(c("Generating validation datasets ... \n"))
  valid.set <- GetTrainTest(valid.train.y, valid.test.y)
  cat(c("Generating evaluation datasets ... \n"))
  eval.set <- GetTrainTest(train.y, test.y)
  return(list(valid=valid.set, eval=eval.set))
}

GetTrainTest <- function(train.years, test.year) {
  # Get interaction and confidence (based on flower abundance) data  
  train.data <- AvailFlwMat(train.years)
  test.data <- AvailFlwMat(test.year)
      
  # Get train datasets
  C.train <- train.data$C.data # interaction frequencies
  A.train.mat <- train.data$abund.mat # flower abundance over time and space  
  like.train <- train.data$like.mat # 
  dislike.train <- train.data$dislike.mat
  abund.train <- train.data$abund.mat       
  flw.occur.train <- train.data$flw.occur 
  A.train.vec <- GetFlwAbundGivenPlants(train.data$flw.data, colnames(C.train), 
                                        FALSE) # avergaged flower abundance
  
  # Remove noisy data 
  remove.idx <- which(A.train.vec==0)
  if(length(remove.idx)>0) {
    #print("Plants that are not in flower survey although it has interactions")
    #print(colnames(C.train)[remove.idx])
    C.train <- C.train[ , -remove.idx]
    A.train.vec <- as.matrix(A.train.vec[-remove.idx, ])
    A.train.mat <- A.train.mat[ , -remove.idx]
    like.train <- like.train[ , -remove.idx]
    dislike.train <- dislike.train[ , -remove.idx]  
    abund.train <- abund.train[ , -remove.idx] 
    flw.occur.train <- flw.occur.train[-remove.idx]
  }

  # Get test datasets
  C.test <- test.data$C.data
  like.test <- test.data$like.mat
  dislike.test <- test.data$dislike.mat
  abund.test <- test.data$abund.mat 
  flw.occur.test <- test.data$flw.occur 
 
  # Reduce train & test matrix with "common" species
  common.test.mats <- GetCommonMatrics(C.train, C.test, like.test, dislike.test, 
                                       abund.test, flw.occur.test)  
  C.train.common <- common.test.mats[[1]]    
  C.test.common <- common.test.mats[[2]]   
  like.test.common <- common.test.mats[[3]]
  dislike.test.common <- common.test.mats[[4]]   
  A.test.common <- common.test.mats[[5]]  
  flw.occur.common <- common.test.mats[[6]]    
  C.test.new <- GetNewInt(C.train.common, C.test.common)    
    
  A.test.vec <- GetFlwAbundGivenPlants(test.data$flw.data, 
                                       colnames(C.test.common), FALSE) 
  
  # Remove noisy data 
  remove.idx <- which(A.test.vec==0)
  if(length(remove.idx)>0) {
    A.test.vec <- as.matrix(A.test.vec[-remove.idx, ])
    C.test.common <- C.test.common[ , -remove.idx]
    like.test.common <- like.test.common[ , -remove.idx]
    dislike.test.common <- dislike.test.common[ ,-remove.idx]        
    A.test.common <- A.test.common[ , -remove.idx]
    abund.test <- abund.test[ , -remove.idx]  
    flw.occur.common <- flw.occur.common[-remove.idx]
    C.test.new <- C.test.new[ , -remove.idx] 
  }  
  
  return(list(Y=list(train=train.years, test=test.year),
              C=list(train=C.train, test=C.test.common), 
              N=list(test=C.test.new),
              Avec=list(train=A.train.vec, test=A.test.vec), 
              Amat=list(train=A.train.mat, test=A.test.common), 
              Like=list(train=like.train, test=like.test.common), 
              Dislike=list(train=dislike.train, test=dislike.test.common),
              Occur=list(train=flw.occur.train, test=flw.occur.common)))
}

# ============================================================================ #
# User-defined datasets
# mydata
# - valid
#   - Y # train/test years
#     - train
#     - test
#   - C # interaction matrix (interaction counts)
#     - train
#     - test
#   - N # indicate matrix for newly appeared interactions
#     - test
#   - A # flower abundance at the test phase
#     - test
#   - R # positive implicit feedback (confidence)
#     - train
#   - D # negative implicit feedback (confidence)
#     - train  
# - eval
#   ... as same as the valid set
# ============================================================================ #
UserDefinedDatasets <- function(test.y) {
  # year configuration
  if(test.y == year.list[1]){
    valid.test.y <- tail(year.list, n=1)
  } else {
    valid.test.y <- test.y - 1
  }  
  # Set the training set for validation
  valid.train.y <- year.list[!year.list== valid.test.y & !year.list== test.y]
  # Set the training set for testing
  train.y <- year.list[!year.list==test.y]
  valid.Y = list(train=valid.train.y, test=valid.test.y)
  
  # simulating datasets
  # For validation,
  valid.C = list(train=matrix(sample(15,100,T),10, dimnames = list(1:10, 1:10)), 
    test=matrix(sample(15,25,T),5, dimnames = list(seq(1,10,2), seq(1,10,2))))
  valid.N = list(test=matrix(sample(c(0,1),25,T),5, dimnames = list(seq(1,10,2), seq(1,10,2))))
  valid.A = list(test=matrix(sample(15,25,T),5, dimnames = list(seq(1,10,2), seq(1,10,2))))
  valid.R = list(train=matrix(sample(15,100,T),10, dimnames = list(1:10, 1:10)))
  valid.D = list(train=matrix(sample(15,100,T),10, dimnames = list(1:10, 1:10)))
  myValid = list(Y=valid.Y, C=valid.C, N=valid.N, A=valid.A, R=valid.R, D=valid.D)
  
  # For evaluation,
  eval.Y = list(train=train.y, test=test.y)
  eval.C = list(train=matrix(sample(15,400,T),20, dimnames = list(1:20, 1:20)), 
    test=matrix(sample(15,100,T),10, dimnames = list(seq(1,20,2), seq(1,20,2))))
  eval.N = list(test=matrix(sample(15,100,T),10, dimnames = list(seq(1,20,2), seq(1,20,2))))
  eval.A = list(test=matrix(sample(15,100,T),10, dimnames = list(seq(1,20,2), seq(1,20,2))))
  eval.R = list(train=matrix(sample(15,400,T),20, dimnames = list(1:20, 1:20)))
  eval.D = list(train=matrix(sample(15,400,T),20, dimnames = list(1:20, 1:20)))
  myEval= list(Y=eval.Y, C=eval.C, N=eval.N, A=eval.A, R=eval.R, D=eval.D)

  mydata = list(valid=myValid, eval=myEval)
  
  return(mydata)
}