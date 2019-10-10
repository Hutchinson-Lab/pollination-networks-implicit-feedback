# ============================================================================ #
# A closed form of matrix factorization (Lin, Kamar, and Horvitz et al., 2014)
# Author: Rebecca A. Hutchinson 
# Data: February 24, 2018
# ============================================================================ #

logistic <- function(x,a,b) { (1 + exp(a*(log(x)+b)))^-1 }

fitUV_C <- function(R,W,lambda,k) {  
  nPolls <- dim(R)[1]
  nPlants <- dim(R)[2]
  P <- (R>0)*1

  nIters <- 20  
  #U <- matrix(runif(nPolls*k),nrow=nPolls,ncol=k)
  #V <- matrix(runif(nPlants*k),nrow=k,ncol=nPlants)
  U <- matrix(0.1,nrow=nPolls,ncol=k)
  V <- matrix(0.1,nrow=k,ncol=nPlants)

  for(n in 1:nIters) {
    oldU <- U
    oldV <- V
    for (i in 1:nPolls) {
      U[i,] <- solve(V%*%diag(W[i,])%*%t(V) + 
                     diag(lambda,k,k))%*%V%*%diag(W[i,])%*%P[i,]
    }
    for (j in 1:nPlants) {
      V[,j] <- solve(t(U)%*%diag(W[,j])%*%U + 
                     diag(lambda,k,k))%*%t(U)%*%diag(W[,j])%*%P[,j]
    }
    #print(paste("Iter ", n, " mean absolute changes in factor U: ", 
    #      mean(abs(oldU-U)), " and factor V: ", mean(abs(oldV-V)), sep=""))
  }
  
  pred_mat <- U%*%V 
  #saveRDS(U, file=paste("RData/U1.RData", sep=""))
  #saveRDS(V, file=paste("RData/V1.RData", sep=""))
  
  pred_mat <- Normalize(pred_mat)
  return(pred_mat)
}