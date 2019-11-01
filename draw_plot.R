# ==================== #
# Draw PR cureves                                                                
# Author: Eugene Seo 
# Data: Sep. 25, 2017
# ==================== #

colors <- rainbow(7)
linetype <- c(2,2,2,2,1,1,1)
plotchar <- c(2,4,5,20,15,17,8)
linewidth <- c(4,4,4,4,6,6,6)
mycex <- 3.5

PRCurve <- function(LB, Pop, Avail, ItemNN, UserNN, MF, IFMF, IFMF2, y, t, v) {
  xrange <- range(Pop[,2], Avail[,2], ItemNN[,2], UserNN[,2], MF[,2], IFMF[,2], IFMF2[,2])
  yrange <- range(Pop[,1], Avail[,1], ItemNN[,1], UserNN[,1], MF[,1], IFMF[,1], IFMF2[,1])
  
  plot(xrange, yrange, type="n", xlab='', ylab='', cex.axis=mycex+1)

  # add lines  
  lines(Pop[,2], Pop[,1], type="b", lwd=linewidth[1], lty=linetype[1], col=colors[1], pch=plotchar[1])
  lines(Avail[,2], Avail[,1], type="b", lwd=linewidth[2], lty=linetype[2], col=colors[2], pch=plotchar[2])  
  lines(ItemNN[,2], ItemNN[,1], type="b", lwd=linewidth[3], lty=linetype[3], col=colors[3], pch=plotchar[3])
  lines(UserNN[,2], UserNN[,1], type="b", lwd=linewidth[4], lty=linetype[4], col=colors[4], pch=plotchar[4])
  lines(MF[,2], MF[,1], type="b", lwd=linewidth[5], lty=linetype[5], col=colors[5], pch=plotchar[5])
  lines(IFMF[,2], IFMF[,1], type="b", lwd=linewidth[6], lty=linetype[6], col=colors[6], pch=plotchar[6])
  lines(IFMF2[,2], IFMF2[,1], type="b", lwd=linewidth[7], lty=linetype[7], col=colors[7], pch=plotchar[7])
  grid(lty=5, lwd=3, col=gray(.9))
  
  # add a title and subtitle
  if(v == 3 | v == 1) {
    atype <- "vec"
  } else if(v == 4 | v == 2) {
    atype <- "mat"
  } else {
    atype <- "sim"
  }  
  title(paste(t, " interactions; A.", atype, sep=""), cex.main=mycex+1.7)
}