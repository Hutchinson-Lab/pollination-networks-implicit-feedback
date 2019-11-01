#==============================================================================#
# Description: Plant-pollinators Interaction Matrix Analysis
# Author: Eugene Seo
# Date: April 30, 2018
#==============================================================================#

LoadMatrix <- function(){
  myData <- LoadData()
  myMat <- BuildPPIMat(myData$int.data, FALSE)
  return(myMat)
}

library(reshape)
BuildPPIMat <- function(data, binary) {
  int.mat <- cast(data, PLTSP_NAME~VISSP_NAME, value = "NO_INT",
                  fun.aggregate = sum, na.rm = TRUE)

  if (dim(int.mat)[1] == 1 & dim(int.mat)[2] == 2) {
    dim.names <- list(int.mat[1,1], names(int.mat)[2])
    int.mat <- matrix(int.mat[1,2], dimnames = dim.names)
  } else {
    int.mat <- as.matrix(int.mat)
  }

  if (binary == TRUE) {
    int.mat <- (int.mat > 0) * 1
  }

  int.mat <- t(int.mat) # format of pollinators x plants
  #cat(c("Matrix size:", dim(int.mat), "\n"))
  #cat(c("Matrix density:", Connectance(int.mat), "\n"))
  return(int.mat)
}

SortMatrix <- function(myMatrix) {
  row.sum = rowSums(myMatrix)
  col.sum = colSums(myMatrix)
  row.order = order(-row.sum)
  col.order = order(-col.sum)
  sortedMat = myMatrix[row.order,]
  sortedMat = sortedMat[,col.order]
  return(sortedMat)
}

FilterMatrixByPerc <- function(sortedMat, perc) {
  total.sum = sum(sortedMat)
  threshold = total.sum * perc
  row.sum = rowSums(sortedMat)
  col.sum = colSums(sortedMat)
  row.cum = cumsum(row.sum)
  col.cum = cumsum(col.sum)
  row.cut = length(which(row.cum < threshold))
  col.cut = length(which(col.cum < threshold))
  subMat <- sortedMat[1:row.cut, 1:col.cut]
  cat(c("Filtered matrix size:", dim(subMat), "\n"))
  return(subMat)
}

FilterMatrixByCounts <- function(sortedMat, counts) {
  row.sum = rowSums(sortedMat)
  row.cut = length(which(row.sum > counts))
  subMat <- sortedMat[1:row.cut, ]
  col.sum = colSums(subMat)  
  remove.idx = which(col.sum == 0)
  subMat <- subMat[, -remove.idx]
  cat(c("Filtered matrix size:", dim(subMat), "\n"))
  return(subMat)
}

# Connectance = density
Connectance <- function(mat) {
  mat <- (mat>0)*1
  return(sum(mat) * 100 / length(mat))
}