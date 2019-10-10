#==============================================================================#
# Description: Manipulate or query EISI PPI data
# Author: Eugene Seo
# Date: February 24, 2018
#==============================================================================#

# ---------------------------------------------------------------------------- #
# --- Load EISI PPI data without NA entries ---------------------------------- #
# ---------------------------------------------------------------------------- #
LoadData <- function() {
  # Load PPI data and flower data from files
  int.data <- read.csv("Data/official_dataset_11-15/SA02601_v2.csv", na.strings = c("","NA"))
  flw.data <- read.csv("Data/official_dataset_11-15/SA02602_v2.csv", na.strings = c("","NA"))

  # Remove data that have NA for plant (PLTSP_NAME) and pollinator (VISSP_NAME)
  # names
  int.data <- RemoveNAInteractionData(int.data)
  flw.data <- RemoveNAFlowerData(flw.data)

  # Filter out data from sub meadows (BD, BH, BS, CNE, CNM, CNT)
  int.data <- FilterOutSubMeadows(int.data)
  flw.data <- FilterOutSubMeadows(flw.data)
  return(list(int.data = int.data, flw.data = flw.data))
}

# Remove interaction data that has NA values for plant and pollinator names
RemoveNAInteractionData <- function(int.data) {
  na <- which(is.na(int.data$PLTSP_NAME) | is.na(int.data$VISSP_NAME))
  if (length(na) > 0) {
    int.data <- int.data[-na, ]
  }
  return(int.data)
}

# Remove flower data that has NA values for species name and abundance features
RemoveNAFlowerData <- function(flw.data) {
  na <- which(is.na(flw.data$PLTSP_NAME) | is.na(flw.data$NO_STALK) |
             is.na(flw.data$NO_FLWS))
  if (length(na) > 0) {
    flw.data <- flw.data[-na, ]
  }
  return(flw.data)
}

# Print the total number of collected data
PrintTotalObservationNum <- function(int.data) {
  print("The number of collected data")
  print(nrow(int.data))
}

# Remove sub meadows which are not relevant or insignificant
FilterOutSubMeadows <- function(input.data) {
  exclusive.meadows <- c("BD", "BH", "BS", "CNE", "CNM", "CNT")
  input.data <- input.data[!input.data$MEADOW%in%exclusive.meadows, ]
  return(input.data)
}

# ---------------------------------------------------------------------------- #
# Get sub interaction and flower data by keywords such as year, meadow, watch
# ---------------------------------------------------------------------------- #
# Get sub interaction & flower data by year directly from original data
# e.g. years = c(11, 12, 13, 14, 15, ...)
GetSubDataByYear <- function(years) {
  myData <- LoadData()
  int.data <- GetSubIntByYear(myData$int.data, years)
  flw.data <- GetSubFlwByYear(myData$flw.data, years)
  return(list(int.data=int.data, flw.data=flw.data))
}

# Get sub interaction data by year from given data
GetSubIntByYear <- function(int.data, years) {
  years <- lapply(years, function(years) paste("20", years, sep=""))
  if (length(years) > 1) {
    sub.int <- int.data[int.data$YEAR%in%years, ]
  } else {
    sub.int <- int.data[int.data$YEAR==years, ]
  }
  return(sub.int )
}

# Get sub flower data by year from given data
GetSubFlwByYear <- function(flw.data, years) {
  years <- lapply(years, function(years) paste("20", years, sep=""))
  if (length(years) > 1) {
    sub.flw <- flw.data[flw.data$YEAR%in%years, ]
  } else {
    sub.flw <- flw.data[flw.data$YEAR==years, ]
  }
  return(sub.flw)
}

# Get sub interaction & flower data by year & meadow directly from origin data
GetSubDataByYearMeadow <- function(y, m) {
  myData <- LoadData()
  int.data <- GetSubIntByYearMeadow(myData$int.data, y, m)
  flw.data <- GetSubFlwByYearMeadow(myData$flw.data, y, m)
  return(list(int.data=int.data, flw.data=flw.data))
}

# Get sub interaction data by year & meadow from given data
GetSubIntByYearMeadow <- function(int.data, y, m) {
  y <- paste("20", y, sep="")
  sub.int.data <- int.data[int.data$YEAR==y & int.data$MEADOW==m, ]
  return(sub.int.data)
}

# Get sub flower data by year & meadow from given data
GetSubFlwByYearMeadow <- function(flw.data, y, m) {
  y <- paste("20", y, sep="")
  sub.flw.data <- flw.data[flw.data$YEAR==y & flw.data$MEADOW==m, ]
  return(sub.flw.data)
}

# Get sub interaction & flower data by year & meadow & watch from original data
GetSubDataByYearMeadowWatch <- function(y,m,w) {
  myData <- LoadData()
  int.data <- GetSubIntByYearMeadowWatch(myData$int.data, y, m, w)
  flw.data <- GetSubFlwByYearMeadowWatch(myData$flw.data, y, m, w)
  return(list(int.data=int.data, flw.data=flw.data))
}

# Get sub interaction data by year & meadow & watch from given data
GetSubIntByYearMeadowWatch <- function(int.data, y, m, w) {
  y <- paste("20", y, sep="")
  sub.int.data <- int.data[int.data$YEAR==y & int.data$MEADOW==m &
                           int.data$WATCH==w, ]
  return(sub.int.data)
}

# Get sub flower data by year & meadow & watch from given data
GetSubFlwByYearMeadowWatch <- function(flw.data, y, m, w) {
  y <- paste("20", y, sep="")
  sub.flw.data <- flw.data[flw.data$YEAR==y & flw.data$MEADOW==m &
                           flw.data$WATCH==w, ]
  return(sub.flw.data)
}

GetCommonMatrics <- function(mat1, mat2, mat3=FALSE, mat4=FALSE, mat5=FALSE, mat6=FALSE) {
  common.pollis <- GetCommon(rownames(mat1), rownames(mat2))
  common.plants <- GetCommon(colnames(mat1), colnames(mat2))
  mat1 <- mat1[rownames(mat1) %in% common.pollis,
               colnames(mat1) %in% common.plants]
  mat2 <- mat2[rownames(mat2) %in% common.pollis,
               colnames(mat2) %in% common.plants]
  if(is.matrix(mat3)) {
    mat3 <- mat3[rownames(mat3) %in% common.pollis,
                 colnames(mat3) %in% common.plants]
    mat4 <- mat4[rownames(mat4) %in% common.pollis,
                 colnames(mat4) %in% common.plants]
    mat5 <- mat5[rownames(mat5) %in% common.pollis,
                 colnames(mat5) %in% common.plants]
    mat6 <- mat6[rownames(mat6) %in% common.plants]
  }
  return(list(mat1, mat2, mat3, mat4, mat5, mat6))
}

GetCommon <- function(data1, data2) {
  return(intersect(data1, data2))
}

# Extract a test set with the consideration of new interactions
GetNewInt <- function(mat1, mat2) {
  mat1 <- (mat1 > 0) * 1
  mat2 <- (mat2 > 0) * 1
  mat3 <- mat2 - mat1
  mat3[which(mat3 != 1)] <- 0
  return(mat3)
}

# Row normalization 0 to 1
Normalize <- function(x) {
  #return ((x - min(x)) / (max(x) - min(x)))
  return (t(apply(x, 1, function(x)(x-min(x))/(max(x)-min(x)))))
  #return (t(apply(x, 1, function(x)(x)/(max(x)))))
}