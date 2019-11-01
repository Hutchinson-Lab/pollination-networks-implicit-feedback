# ============================================================================ #
# Extract various formats of flower availability
#  1. flw.data: The total number of flowers that were available 
#  2. like.mat: The number of times a plant and a pollinator were interacted
#  3. dislike.mat: The number of times a plant and a pollinator were together 
#                  but not related      
#  4. chance.mat: The number of a plant and a pollinator were together    
#  5. flw.occur: The number of a plant appeared in all surveys                     
# Author: Eugene Seo 
# Data: February 24, 2018
# ============================================================================ # 

AvailFlwMat <- function(years){
  mydata <- GetSubDataByYear(years)      
  C.data <- BuildPPIMat(mydata$int.data, FALSE) # interaction frequencites
  abund.mat <- matrix(0, nrow=dim(C.data)[1], ncol=dim(C.data)[2])
  chance.mat <- matrix(0, nrow=dim(C.data)[1], ncol=dim(C.data)[2])
  like.mat <- matrix(0, nrow=dim(C.data)[1], ncol=dim(C.data)[2])
  dislike.mat <- matrix(0, nrow=dim(C.data)[1], ncol=dim(C.data)[2])
  flw.occur <- array(0, c(ncol(C.data)))
  plant.occur <- array(0, c(ncol(C.data))) #including not interacting ones
  poll.occur <- array(0, c(nrow(C.data)))  

  # define row and column names 
  poll.names = rownames(dislike.mat) = rownames(poll.occur) = rownames(C.data)
  rownames(abund.mat) = rownames(chance.mat) = rownames(like.mat) = poll.names
  plant.names = rownames(flw.occur) = rownames(plant.occur) = colnames(C.data)  
  colnames(abund.mat) = colnames(chance.mat) = plant.names
  colnames(like.mat) = colnames(dislike.mat) = plant.names
  
  count <- 0
  for(y in years) {
    cat(c("\tBuilding confidence matrix for", y, "year ...\n"))
    year.data <- GetSubDataByYear(y)    
    meadows <- sort(unique(year.data$int.data$MEADOW))
    for(m in meadows) {
      meadow.data <- year.data$int.data[year.data$int.data$MEADOW==m, ]
      watches <- sort(unique(meadow.data$WATCH))
      for(w in watches) {
        # Get sub-interactions that occured at the specific watch & meadow
        sub.data <- GetSubDataByYearMeadowWatch(y, m, w)
        sub.mat <- BuildPPIMat(sub.data$int.data, FALSE)
        P <- (sub.mat>0)*1
        flw.data <- sub.data$flw.data[order(sub.data$flw.data$PLTSP_NAME), ]
        
        # Get common (interacted) plants' flower abundnace
        plants <- intersect(colnames(C.data), sub.data$flw.data$PLTSP_NAME)
        avail <- GetFlwAbundGivenPlants(flw.data, plants, FALSE)        
        pollis <- unique(sub.data$int.data$VISSP_NAME)

        # Remove noisy data (plants without flowers in the interaction matrix)
        remove.idx <- which(!colnames(sub.mat) %in% plants)
        if(length(remove.idx) > 0) {
          sub.mat <- sub.mat[,-remove.idx]
          P <- P[,-remove.idx]
        }
        
        if(is.matrix(sub.mat)) {
          row.idx <- rownames(C.data) %in% pollis
          col.idx <- colnames(C.data) %in% plants 
          # The flower abundance at the meadow & watch (MW)
          avail.mat = t(matrix(rep(avail,length(pollis)), ncol=length(pollis)))
          # Summed MW-level flower abundance over all meadows and watches (MWs)
          abund.mat[row.idx, col.idx] <- abund.mat[row.idx, col.idx] + avail.mat
          # Summed co-occurrence number of flowers over all MWs for pollinators
          chance.mat[row.idx, col.idx] <- chance.mat[row.idx, col.idx] + 1
          # number of appearance of flowers over MWs
          flw.occur[col.idx] <- flw.occur[col.idx] + 1
          col.idx <- colnames(C.data) %in% colnames(sub.mat)
          # Confidence of likes (Positive implicit feedback)
          # = number of MWs where interactions occured
          if((length(which(row.idx))*length(which(col.idx))) == length(P)){
            like.mat[row.idx, col.idx] <- like.mat[row.idx, col.idx] + P
          } else {
            print("ERROR 1")
            cat(years, y, m, w)
          }
        }        
        mw.poll.names = unique(sub.data$int.data$VISSP_NAME)
        mw.plant.names = unique(sub.data$flw.data$PLTSP_NAME)
        poll.idx = which(poll.names %in% mw.poll.names)
        plant.idx = which(plant.names %in% mw.plant.names)
        poll.occur[poll.idx] = poll.occur[poll.idx] + 1
        plant.occur[plant.idx] = plant.occur[plant.idx] + 1
        count <- count + 1 	 # total number of meadow-watches (MWs)
      }
    }
  }
  # Condifence of dislikes (Negative implicit feedback)
  # = number of MWs where no interactions occured in the same contexts
  if(length(which(dislike.mat<0))>0){
    print("ERROR 2")
  }
  dislike.mat <- chance.mat - like.mat 
  
  #freq.mat <- like.mat / chance.mat

  return(list(C.data=C.data, flw.data=mydata$flw.data, like.mat=like.mat, 
              dislike.mat=dislike.mat, abund.mat=abund.mat, 
              chance.mat=chance.mat, flw.occur=flw.occur, 
              plant.occur=plant.occur, poll.occur=poll.occur, count=count))
}

# Extract flower abundance from given plant species names
GetFlwAbundGivenPlants <- function(flw.data, plants.names, logsum=FALSE) {
  plant.num <- length(plants.names)
  dim.names <- list(plants.names, c("Tot.Flw"))
  flw.abund <- matrix(0, nrow=plant.num, ncol=1, dimnames=dim.names)
  for (pl in 1:plant.num) {
    idx <- which(flw.data$PLTSP_NAME==plants.names[pl])
    flw.abund.list = flw.data$NO_STALK[idx] * flw.data$NO_FLWS[idx]
    if (length(idx)>0) {
      if(logsum) {
        flw.abund[pl] <- sum(log(flw.abund.list), na.rm=T)
      }
      else {
        flw.abund[pl] <- sum(flw.abund.list, na.rm=T)
      }
    }
  }
  return(flw.abund)
}