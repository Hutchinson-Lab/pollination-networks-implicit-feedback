# ============================================================================ #
# Main Functions
# Usages:
#  1. Generate data of a specific year: e.g. GenerateDataset(11)
#  2. Tune parameteres: e.g. Validation(11)
#  3. Learn and evaluate models: e.g. Evaluation(11)
# Author: Eugene Seo
# Data: February 24, 2018
# ============================================================================ #
remove(list=ls())
source("data_processing.R")
source("data_generation.R")
source("cross_validation.R")
source("matrix_analysis.R")
source("baseline.R")
source("draw_plot.R")


# ---------------------------------------------------------------------------- #
# Experiment factors
# 1. Version: Different format of utilizing availability info.
#    1:AV-AV, 2:AM-AM, 3:OV-OV, 4:OM-OM, 5:AM-AV, 6:OM-AV, 7:OM-AM
#      - A: Flower availability
#      - O: Number of opportunity to interact with flowers over surveys
#      - V: Vector format,
#      - M: Matrix format (personalized format)
# 2. Target: Type of target links to be predicted
#  - All: all test interactions
#  - New: test interactions appearing only in test set but not in training set
# 3. Year: Survey years
# 4. Model: List of models to be used for experiments.
# ---------------------------------------------------------------------------- #
version.list <- c(3, 4)
target.list <- c("All", "New")
year.list <- c(11,12,13,14,15)
model.list <- c("LB", "Pop", "Avail", "ItemNN", "UserNN", "MF", "IFMF", "IFMF2")


# 1. Generate training and test sets
GenerateDataset <- function(test.y) {
  cat("Generating train/test sets for the test year", test.y, " ...\n")
  mydata <- GetValidEval(test.y)
  # DATA SAVE
  if (!dir.exists('RData')) {
    dir.create('RData')
  }
  saveRDS(mydata, file=paste("RData/mydata_", test.y, ".RData", sep=""))
}

# 2. Cross valdiation on validation sets
Validation <- function(test.y) {
  if (!dir.exists('Param')) {
    dir.create('Param')
  }
  cat("Cross Valdiation...", test.y, "\n")
  # DATA LOAD
  filename <- paste("RData/mydata_", test.y, ".RData", sep="")
  if(file.exists(filename)){
    mydata <- readRDS(paste("RData/mydata_", test.y, ".RData", sep=""))
  } else {
    mydata = SimulateDataset(test.y)
    version.list = c("Sim")
  }

  for(ver in version.list) {
    for(target in target.list) {
      df <- CrossValidation(mydata$valid, test.y, ver, target)
      # PARAM SAVE
      filename <- paste("Param/Param_", ver, "_", target, ".csv", sep="")
      if (!file.exists(filename)) {
        write.table(df, filename, sep=",", col.names=NA, append=TRUE)
      } else {
        write.table(df, filename, sep=",", col.names=FALSE, append=TRUE)
      }
    }
  }
}

# 3. Evaluation models
Evaluation <- function(test.y) {
  if (!dir.exists('Results')) {
    dir.create('Results')
  }
  
  cat("Testing... test.y", test.y, "\n")
  # DATA LOAD
  filename <- paste("RData/mydata_", test.y, ".RData", sep="")
  if(file.exists(filename)){
    mydata <- readRDS(filename)
  } else {
    mydata = SimulateDataset(test.y) 
    version.list = c("Sim")
  }
  
  attach(mtcars)
  x11(width=30, height=20, pointsize=20)
  png(paste("Results/PR_34_", test.y, ".png", sep=""),
            width=20, height=20, units="in", res=300, pointsize=15)
  par(mfrow=c(length(target.list),length(version.list)))

  for(target in target.list) {
    df <- data.frame(matrix(NA, ncol=16)); idx <- 1
    for(ver in version.list) {
       cat("\n[test.y:", test.y, "version:", ver, "target:", target, "]\n")
       # PARAM LOAD
       params <- LoadTunedParam(ver, target)
       opt.params <- GetTunedParamsByYear(params, test.y)
       if (dim(opt.params)[1] > 1) {
         print("Please remove previous param files")
         stop()
       }
       res <- EvaluateTestset(mydata$eval, opt.params, ver, target, test.y,
                              plot.type)
       # RESULT SAVE
       saveRDS(res, file=paste("Results/PR_", test.y, "_", ver, "_", target,
                               ".RData", sep=""))
       df <- WriteResult(df, idx, res, mydata$eval, model.list, ver, target);
       idx <- idx + 1
    }
    colnames(df) <- c("version", "train years", "train size", "test year",
                      "test size", "test int. type", "int. num", model.list)
    # RESULT SAVE
    write.table(df, paste("Results/MPR_", target, ".csv", sep=""), col.names=NA,
                sep=",", append=TRUE)
  }
  dev.off()
}

EvaluateTestset <- function(mydata, opt.params, ver, target, test.y, plot.type) {
  C <- mydata$C
  N <- mydata$N
  P <- (C$train>0)*1
  nPolls <- nrow(C$train)
  nPlant <- ncol(C$train)

  if(is.null(mydata$A)) {
    Avec <- mydata$Avec
    Amat <- mydata$Amat
    Like <- mydata$Like
    Dislike <- mydata$Dislike
    Occur <- mydata$Occur 
    
    if(ver == 1) { # AV-AV
      R <- C$train
      D <- t(matrix(rep(Avec$train, nPolls), ncol=nPolls))
      Atest <- t(matrix(rep(Avec$test,nrow(C$test)),ncol=nrow(C$test)))
    } else if(ver == 2) { # AM-AM
      R <- C$train
      D <- Amat$train
      Atest <- Amat$test
    } else if(ver == 3) { # OV-OV
      R <- C$train
      D <- t(matrix(rep(Occur$train, nPolls), ncol=nPolls))
      Atest <- t(matrix(rep(Occur$test,nrow(C$test)),ncol=nrow(C$test)))
    } else if(ver == 4) { # OM-OM
      R <- Like$train
      D <- Dislike$train
      Atest <- Like$test + Dislike$test
    } else if(ver == 5) { # AM-AV
      R <- C$train
      D <- Amat$train
      Atest <- t(matrix(rep(Avec$test,nrow(C$test)),ncol=nrow(C$test)))
    } else if(ver == 6) { # OM-AV
      R <- Like$train
      D <- Dislike$train
      Atest <- t(matrix(rep(Avec$test,nrow(C$test)),ncol=nrow(C$test)))
    } else if(ver == 7) { # OM-AM
      R <- Like$train
      D <- Dislike$train
      Atest <- Amat$test
    }
  } else {
    R <- mydata$R$train
    D <- mydata$D$train
    Atest <- mydata$A$test
  }

  # LB
  cat("LB", "\n")
  LB <- EvaluationAll(C$test, C$test, N$test, target)

  # Popularity
  cat("Popularity", "\n")
  int.mat <- Baseline(C$train, "popularity")
  pred.mat <- Normalize(int.mat)
  Pop <- EvaluateCommon(C$train, C$test, N$test, Atest, pred.mat, target)

  # Availability
  cat("Availability", "\n")
  Avail <- EvaluationAll(C$test, Atest, N$test, target)

  # Neighborhood
  cat("ItemNN... params", opt.params$ITEM.k, "\n")
  int.mat <- Baseline(C$train, "ItemNN", opt.params$ITEM.k)
  pred.mat <- Normalize(int.mat)
  ItemNN <- EvaluateCommon(C$train, C$test, N$test, Atest, pred.mat, target)

  cat("UserNN... params", opt.params$USER.k, "\n")
  int.mat <- Baseline(C$train, "UserNN", opt.params$USER.k)
  pred.mat <- Normalize(int.mat)
  UserNN <- EvaluateCommon(C$train, C$test, N$test, Atest, pred.mat, target)

  # MF
  cat("MF... params", opt.params$MF.lam, opt.params$MF.k, "\n")
  W0 <- matrix(1, nrow=nPolls, ncol=nPlant)
  pred.mat <- fitUV_C(C$train, W0, opt.params$MF.lam, opt.params$MF.k)
  MF <- EvaluateCommon(C$train, C$test, N$test, Atest, pred.mat, target)

  # IFMF
  cat("IFMF... params", opt.params$IFMF.a, opt.params$IFMF.lam,
                        opt.params$IFMF.k, "\n")
  W1 <- 1 + opt.params$IFMF.a * R
  W1 <- W1/max(W1)
  pred.mat <- fitUV_C(C$train, W1, opt.params$IFMF.lam, opt.params$IFMF.k)
  IFMF <- EvaluateCommon(C$train, C$test, N$test, Atest, pred.mat, target)

  # IFMF2
  cat("IFMF2...params", opt.params$IFMF2.a1, opt.params$IFMF2.b1,
                        opt.params$IFMF2.a2, opt.params$IFMF2.b2,
                        opt.params$IFMF2.lam, opt.params$IFMF2.k, "\n")
  W2 <- W0
  W2[P==1] <- logistic(R[P==1], opt.params$IFMF2.a1, opt.params$IFMF2.b1)
  W2[P==0] <- logistic(D[P==0], opt.params$IFMF2.a2, opt.params$IFMF2.b2)
  pred.mat <- fitUV_C(C$train, W2, 0.01, opt.params$IFMF2.k)
  IFMF2 <- EvaluateCommon(C$train, C$test, N$test, Atest, pred.mat, target)

  PRCurve(LB$PR, Pop$PR, Avail$PR, ItemNN$PR, UserNN$PR, MF$PR, IFMF$PR,
          IFMF2$PR, mydata$Y$test, target, ver)

  return(list(LB=LB, Pop=Pop, Avail=Avail, ItemNN=ItemNN, UserNN=UserNN, MF=MF,
              IFMF=IFMF, IFMF2=IFMF2))
}

WriteResult <- function(df, idx, res, mydata, model.list, v, target) {
  pos.int.num <- length(which(mydata$C$test>0))
  new.int.num <- length(which(mydata$N$test>0))
  tab <- 7
  if(target == "All") {
    df[idx,1] <- v
    df[idx,2] <- paste(mydata$Y$train, sep="", collapse="")
    df[idx,3] <- paste(dim(mydata$C$train)[1], "x", dim(mydata$C$train)[2],
                       sep="", collapse="")
    df[idx,4] <- mydata$Y$test
    df[idx,5] <- paste(dim(mydata$C$test)[1], "x", dim(mydata$C$test)[2],
                       sep="", collapse="")
    df[idx,6] <- "All"
    df[idx,7] <- pos.int.num
    for(j in 1:length(model.list)) {
      df[idx,(j+tab)] <- res[[j]]$MPR
    }
    j <- j + tab + 1
    df[idx,j] <- min(df[idx, 9:15])
  } else if(target == "New") {
    df[idx,1] <- v
    df[idx,2] <- paste(mydata$Y$train, sep="", collapse="")
    df[idx,3] <- paste(dim(mydata$C$train)[1], "x", dim(mydata$C$train)[2],
                       sep="", collapse="")
    df[idx,4] <- mydata$Y$test
    df[idx,5] <- paste(dim(mydata$C$test)[1], "x", dim(mydata$C$test)[2],
                       sep="", collapse="")
    df[idx,6] <- "New"
    df[idx,7] <- new.int.num
    for(j in 1:length(model.list)) {
      df[idx,(j+tab)] <- res[[j]]$MPR
    }
    j <- j + tab + 1
    df[idx, j] <- min(df[idx, 9:15])
  }
  return(df)
}

EvaluateAllYears <- function() {
  for(test.year in c(11,12,13,14,15)) {
    cat(c("test.year", test.year, "\n"))
    GenerateDataset(test.year)
    Validation(test.year)
    Evaluation(test.year)
  }
}

EvaluateOneYear <- function(test.year) {
  if (test.year %in% c(11,12,13,14,15)) {
    cat(c("test.year", test.year, "\n"))
    GenerateDataset(test.year)
    Validation(test.year)
    Evaluation(test.year)
  }
}

SimulateOneYear <- function(test.year) {
  if (test.year %in% c(11,12,13,14,15)) {
    cat(c("test.year", test.year, "\n"))
    #SimulateDataset(test.year) # TODO
    Validation(test.year)
    Evaluation(test.year)
  }
}