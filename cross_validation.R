# ============================================================================ #
# Tune parameters (k, lambda, alpha, a1, a2, b1, b2) on a validation set
# Author: Eugene Seo 
# Data: February 24, 2018
# ============================================================================ #
source("mf.R")
source("evaluation_mertics.R")

# Tune Parameters
k.list <- c(2, 5, 10, 15)
lambda.list <- c(0.01, 0.1, 1)
a.list <- c(1, 10, 20)

CrossValidation <- function(mydata, test.y, ver, target) {  
  cat("\n[ version", ver, "target", target, "]\n")      
  df <- data.frame()
  idx <- 1
  res <- TunningParameters(mydata, ver, target, test.y)
  df[idx,1] <- paste(mydata$Y$train, collapse="")
  df[idx,2] <- mydata$Y$test
  df[idx,3] <- test.y
  df[idx,4] <- res$MF$k
  df[idx,5] <- res$MF$lambda
  df[idx,6] <- res$MF$MPR
  df[idx,7] <- res$IFMF$k
  df[idx,8] <- res$IFMF$alpha
  df[idx,9] <- res$IFMF$lambda
  df[idx,10] <- res$IFMF$MPR
  df[idx,11] <- res$IFMF2$a1
  df[idx,12] <- res$IFMF2$b1
  df[idx,13] <- res$IFMF2$a2
  df[idx,14] <- res$IFMF2$b2
  df[idx,15] <- res$IFMF2$k
  df[idx,16] <- res$IFMF2$lamda
  df[idx,17] <- res$IFMF2$MPR
  names(df) <- c( "train.y", "valid.y", "test.y", 
                  "MF.k", "MF.lam", "MF.MPR", 
                  "IFMF.k", "IFMF.a", "IFMF.lam", "IFMF.MPR",
                  "IFMF2.a1", "IFMF2.b1", "IFMF2.a2", "IFMF2.b2", 
                  "IFMF2.k", "IFMF2.lam", "IFMF2.MPR")
  return(df)
}

TunningParameters <- function(mydata, ver, target, test.y) {
  cat("MF")
  MF <- MF_Tuning(mydata, k.list, lambda.list, ver, target)
  cat("IFMF")
  IFMF <- IFMF_Tuning(mydata, k.list, lambda.list, a.list, ver, target)
  cat("IFMF2")
  IFMF2 <- IFMF2_Tuning(mydata, k.list, lambda.list, ver, target)
  return(list(MF=MF, IFMF=IFMF, IFMF2=IFMF2))
}

# Tunning lambda
MF_Tuning <- function(mydata, k.list, lambda.list, ver, target) {
  C <- mydata$C
  N <- mydata$N
  Avec <- mydata$Avec
  Amat <- mydata$Amat
  Like <- mydata$Like
  Dislike <- mydata$Dislike
  Occur <- mydata$Occur
  
  if(ver == 1) { # AV-AV
    Avalid <- t(matrix(rep(Avec$test,nrow(C$test)),ncol=nrow(C$test)))
  } else if(ver == 2) { # AM-AM
    Avalid <- Amat$test
  } else if(ver == 3) { # OV-OV
    Avalid <- t(matrix(rep(Occur$test,nrow(C$test)),ncol=nrow(C$test)))
  } else if(ver == 4) { # OM-OM 
    Avalid <- Like$test + Dislike$test
  } else if(ver == 5) { # AM-AV
    Avalid <- t(matrix(rep(Avec$test,nrow(C$test)),ncol=nrow(C$test)))
  } else if(ver == 6) { # OM-AV 
    Avalid <- t(matrix(rep(Avec$test,nrow(C$test)),ncol=nrow(C$test)))
  } else if(ver == 7) { # OM-AM 
    Avalid <- Amat$test
  }

  W0 <- matrix(1, nrow=nrow(C$train), ncol=ncol(C$train))
 
  # k tuning
  MPR.list <- array(0, c(length(k.list)))
  for(i in 1:length(k.list)){
    cat(" k", i)
    k <- k.list[i]
    pred.mat <- fitUV_C(C$train, W0, 0.02, k)
    MPR.list[i] <- ValidateCommon(C$train, C$test, N$test, Avalid, pred.mat, 
                                  target)
  }
  k.opt <- k.list[which.min(MPR.list)]

  # lambda tuning
  MPR.list <- array(0, c(length(lambda.list)))
  for(i in 1:length(lambda.list)){
    cat(" lambda", i)
    lambda <- lambda.list[i]
    pred.mat <- fitUV_C(C$train, W0, lambda, k.opt)
    MPR.list[i] <- ValidateCommon(C$train, C$test, N$test, Avalid, pred.mat, 
                                  target)
  }
  lambda.opt <- lambda.list[which.min(MPR.list)]

  return(list(k=k.opt, lambda=lambda.opt, MPR=min(MPR.list)))
}

IFMF_Tuning <- function(mydata, k.list, lambda.list, a.list, ver, target) {
  C <- mydata$C
  N <- mydata$N
  Avec <- mydata$Avec 
  Amat <- mydata$Amat
  Like <- mydata$Like
  Dislike <- mydata$Dislike
  Occur <- mydata$Occur
  P <- (C$train>0)*1

  if(ver == 1) { # AV-AV
    R <- C$train
    Avalid <- t(matrix(rep(Avec$test,nrow(C$test)),ncol=nrow(C$test)))
  } else if(ver == 2) { # AM-AM
    R <- C$train
    Avalid <- Amat$test
  } else if(ver == 3) { # OV-OV
    R <- Like$train
    Avalid <- t(matrix(rep(Occur$test,nrow(C$test)),ncol=nrow(C$test)))
  } else if(ver == 4) { # OM-OM 
    R <- Like$train
    Avalid <- Like$test + Dislike$test
  } else if(ver == 5) { # AM-AV
    R <- C$train
    Avalid <- t(matrix(rep(Avec$test,nrow(C$test)),ncol=nrow(C$test)))
  } else if(ver == 6) { # OM-AV 
    R <- Like$train
    Avalid <- t(matrix(rep(Avec$test,nrow(C$test)),ncol=nrow(C$test)))
  } else if(ver == 7) { # OM-AM 
    R <- Like$train
    Avalid <- Amat$test
  }

  # alpha tuning
  MPR.list <- array(0, c(length(a.list)))
  for(i in 1:length(a.list)){
    cat(" a", i)
    alpha <- a.list[i]
    W1 <- 1 + alpha * R
    W1 <- W1 / max(W1)
    pred.mat <- fitUV_C(C$train, W1, 0.02, 2)
    MPR.list[i] <- ValidateCommon(C$train, C$test, N$test, Avalid, pred.mat, 
                                  target)
  }
  alpha.opt <- a.list[which.min(MPR.list)]

  W1 <- 1 + alpha.opt * R
  W1 <- W1 / max(W1)

  # k tuning
  MPR.list <- array(0, c(length(k.list)))
  for(i in 1:length(k.list)){
    cat(" k", i)
    k <- k.list[i]
    pred.mat <- fitUV_C(C$train, W1, 0.02, k)
    MPR.list[i] <- ValidateCommon(C$train, C$test, N$test, Avalid, pred.mat, 
                                  target)
  }
  k.opt <- k.list[which.min(MPR.list)]

  # lambda tuning
  MPR.list <- array(0, c(length(lambda.list)))
  for(i in 1:length(lambda.list)){
    cat(" lambda", i)
    lambda <- lambda.list[i]
    pred.mat <- fitUV_C(C$train, W1, lambda, k.opt)
    MPR.list[i] <- ValidateCommon(C$train, C$test, N$test, Avalid, pred.mat, 
                                  target)
  }
  lambda.opt <- lambda.list[which.min(MPR.list)]

  return(list(k=k.opt, alpha=alpha.opt, lambda=lambda.opt, MPR=min(MPR.list)))
}

# Tunning a1, b1, a2, b2, lambda
IFMF2_Tuning <- function(mydata, k.list, lambda.list, ver, target) {
  C <- mydata$C
  N <- mydata$N
  Avec <- mydata$Avec 
  Amat <- mydata$Amat
  Like <- mydata$Like
  Dislike <- mydata$Dislike
  Occur <- mydata$Occur
  P <- (C$train>0)*1
  nPolls <- nrow(C$train)
  nPlant <- ncol(C$train)

  if(ver == 1) { # AV-AV
    R <- C$train
    D <- t(matrix(rep(Avec$train, nPolls), ncol=nPolls))
    Avalid <- t(matrix(rep(Avec$test,nrow(C$test)),ncol=nrow(C$test)))
  } else if(ver == 2) { # AM-AM
    R <- C$train
    D <- Amat$train
    Avalid <- Amat$test
  } else if(ver == 3) { # OV-OV
    R <- Like$train
    D <- t(matrix(rep(Occur$train, nPolls), ncol=nPolls)) 
    Avalid <- t(matrix(rep(Occur$test,nrow(C$test)),ncol=nrow(C$test)))
  } else if(ver == 4) { # OM-OM 
    R <- Like$train
    D <- Dislike$train
    Avalid <- Like$test + Dislike$test
  } else if(ver == 5) { # AM-AV
    R <- C$train
    D <- Amat$train
    Avalid <- t(matrix(rep(Avec$test,nrow(C$test)),ncol=nrow(C$test)))
  } else if(ver == 6) { # OM-AV 
    R <- Like$train
    D <- Dislike$train
    Avalid <- t(matrix(rep(Avec$test,nrow(C$test)),ncol=nrow(C$test)))
  } else if(ver == 7) { # OM-AM 
    R <- Like$train
    D <- Dislike$train
    Avalid <- Amat$test
  }
  
  W <- matrix(, nrow=nPolls, ncol=nPlant)
  W[P==1] <- R[P==1] + 1
  W[P==0] <- D[P==0] + 1

  # logistic f parameters (a1, b1, a2, b2) tuning
  x <- seq(max(C$train[P==1]))
  c_param <- GetParamList(x)
  x <- seq(max(D[P==0]))
  d_param <- GetParamList(x)

  total <- length(c_param[[1]]) * length(c_param[[2]])
  a2 <- d_param[[1]][1]
  b2 <- d_param[[2]][1]    
  param_tuning <- data.frame()
  idx <- 1
  for(a1 in c_param[[1]]) {
    for(b1 in c_param[[2]]) { 
      cat(" a1", a1, "b1", b1, "idx", idx, "/", total, "\n")
      W2 <- W
      W2[P==1] <- logistic(R[P==1], a1, b1)
      W2[P==0] <- logistic(D[P==0], a2, b2)
      param_tuning[idx,1] <- a1
      param_tuning[idx,2] <- b1
      pred.mat <- fitUV_C(C$train, W2, 0.02, 2)
      param_tuning[idx,3] <- ValidateCommon(C$train, C$test, N$test, Avalid, 
                                            pred.mat, target)
      idx <- idx + 1
    }
  }
  opt.param <- param_tuning[which.min(param_tuning[,3]),]
  opt.a1 <- opt.param[[1]]
  opt.b1 <- opt.param[[2]]

  total <- length(d_param[[1]]) * length(d_param[[2]])
  param_tuning <- data.frame()
  idx <- 1
  for(a2 in d_param[[1]]) {
    for(b2 in d_param[[2]]) { 
      cat(" a2", a2, "b2", b2, "idx", idx, "/", total, "\n")
      W2 <- W
      W2[P==1] <- logistic(R[P==1], opt.a1, opt.b1)
      W2[P==0] <- logistic(D[P==0], a2, b2)
      param_tuning[idx,1] <- a2
      param_tuning[idx,2] <- b2
      pred.mat <- fitUV_C(C$train, W2, 0.02, 2)
      param_tuning[idx,3] <- ValidateCommon(C$train, C$test, N$test, Avalid, 
                                            pred.mat, target)
      idx <- idx + 1
    }
  }
  opt.param <- param_tuning[which.min(param_tuning[,3]),]
  opt.a2 <- opt.param[[1]]
  opt.b2 <- opt.param[[2]]
  W2 <- W
  W2[P==1] <- logistic(R[P==1], opt.a1, opt.b1)
  W2[P==0] <- logistic(D[P==0], opt.a2, opt.b2)

  # k tuning
  MPR.list <- array(0, c(length(k.list)))
  for(i in 1:length(k.list)){
    cat(" k", i)
    k <- k.list[i]
    pred.mat <- fitUV_C(C$train, W2, 0.02, k)
    MPR.list[i] <- ValidateCommon(C$train, C$test, N$test, Avalid, pred.mat, 
                                  target)
  }
  k.opt <- k.list[which.min(MPR.list)]

  # lambda tuning
  MPR.list <- array(0, c(length(lambda.list)))
  for(i in 1:length(lambda.list)){
    cat(" lambda", i)
    lambda <- lambda.list[i]
    pred.mat <- fitUV_C(C$train, W2, lambda, k.opt)
    MPR.list[i] <- ValidateCommon(C$train, C$test, N$test, Avalid, pred.mat, 
                                  target)
  }
  lambda.opt <- lambda.list[which.min(MPR.list)] 
  
  return(list(a1=opt.a1, b1=opt.b1, a2=opt.a2, b2=opt.b2, k=k.opt, 
              lamda=lambda.opt, MPR=min(MPR.list)))
}

GetParamList <- function(x) {
  for(beta in -10:-1) {
    if(max(logistic(max(x),-1,beta)) > 0.9) {
      break
    }
  }
  beta.list <- seq(beta, -1, by=0.5)
  for(alpha in -10:-1) {
    if(median(logistic(x, alpha, beta)) < 0.6 | 
       max(logistic(max(x), alpha, beta)) < 0.9) {
      break
    }
  }
  alpha.list <- seq(alpha, -1, by=0.5)
  return(list(a1=alpha.list, b1=beta.list))
}

LoadTunedParam <- function(v, t) {
  filename <- paste("Param/Param_", v, "_", t, ".csv", sep="")
  param <- read.csv(filename)
  return(param)
}

GetTunedParamsByYear <- function(params, test.y) {
  opt <- params[which(params["test.y"]==test.y), ]
  return(opt)
}