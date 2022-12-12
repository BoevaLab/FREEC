#!/usr/bin/env Rscript

# One way to run this script is:
# cat makeGraph.R | R --slave --args <*_ratio.txt> [<*_BAF.txt>]
# Ploidy value will be inferred from the ratio file


args <- commandArgs()

BAFfileInd = 0;
ratioFileInd = 0;

#find which argument is Ratio.txt and which BAF.txt:
for (i in c(1:length(args))) {
  if (length(grep("ratio.txt", args[i]))) {
    ratioFileInd = i;
  }
  if (length(grep("BAF", args[i]))) {
    BAFfileInd = i;
  }
}

#------------------------------------------------------

#plot .png for the _ratio.txt file:

if (ratioFileInd) {
  
  #read the file and get ploidy value:
  
  ratio <-read.table(args[ratioFileInd], header=TRUE);
  ratio<-data.frame(ratio)
  ploidy = median (ratio$CopyNumber[which(ratio$MedianRatio>0.8 & ratio$MedianRatio<1.2)], na.rm = T)
  cat (c("INFO: Selected ploidy:", ploidy, "\n"))
  
  #------------------------------------------------------
  
  #Plotting in the log scale:
  offset = 0.01
  
  png(filename = paste(args[ratioFileInd],".log2.png",sep = ""), width = 1180, height = 1180,
      units = "px", pointsize = 20, bg = "white", res = NA)
  plot(1:10)
  op <- par(mfrow = c(5,5))
  
  for (i in c(1:22,'X','Y')) {
    tt <- which(ratio$Chromosome==i)
    if (length(tt)>0) {
      plot(ratio$Start[tt],log2(ratio$Ratio[tt]+offset),xlab = paste ("position, chr",i),ylab = "normalized copy number profile (log2)",pch = ".",col = colors()[88])
      tt <- which(ratio$Chromosome==i  & ratio$CopyNumber>ploidy )
      points(ratio$Start[tt],log2(ratio$Ratio[tt]+offset),pch = ".",col = colors()[136])
      
      
      tt <- which(ratio$Chromosome==i  & ratio$CopyNumber<ploidy & ratio$CopyNumber!= -1)
      points(ratio$Start[tt],log2(ratio$Ratio[tt]+offset),pch = ".",col = colors()[461])
      tt <- which(ratio$Chromosome==i)
      
      #UNCOMMENT HERE TO SEE THE PREDICTED COPY NUMBER LEVEL:
      #points(ratio$Start[tt],log2(ratio$CopyNumber[tt]/ploidy+offset), pch = ".", col = colors()[24],cex=4)
      
    }
    tt <- which(ratio$Chromosome==i)
    
    #UNCOMMENT HERE TO SEE THE EVALUATED MEDIAN LEVEL PER SEGMENT:
    #points(ratio$Start[tt],log2(ratio$MedianRatio[tt]+offset), pch = ".", col = colors()[463],cex=4)
    
  }
  dev.off()
  
  #------------------------------------------------------
  #Plotting in raw ratio values:
  png(filename = paste(args[ratioFileInd],".png",sep = ""), width = 1180, height = 1180,
      units = "px", pointsize = 20, bg = "white", res = NA)
  plot(1:10)
  op <- par(mfrow = c(5,5))
  
  #replace high values of ratio with value "maxLevelToPlot":
  maxLevelToPlot <- 3
  ratio$Ratio[ratio$Ratio>maxLevelToPlot]=maxLevelToPlot
 
  for (i in c(1:22,'X','Y')) {
    tt <- which(ratio$Chromosome==i)
    if (length(tt)>0) {
      plot(ratio$Start[tt],ratio$Ratio[tt]*ploidy,ylim = c(0,maxLevelToPlot*ploidy),xlab = paste ("position, chr",i),ylab = "normalized copy number profile",pch = ".",col = colors()[88])
      tt <- which(ratio$Chromosome==i  & ratio$CopyNumber>ploidy )
      points(ratio$Start[tt],ratio$Ratio[tt]*ploidy,pch = ".",col = colors()[136])
      
      tt <- which(ratio$Chromosome==i  & ratio$Ratio==maxLevelToPlot & ratio$CopyNumber>ploidy)	
      points(ratio$Start[tt],ratio$Ratio[tt]*ploidy,pch = ".",col = colors()[136],cex=4)
      
      tt <- which(ratio$Chromosome==i  & ratio$CopyNumber<ploidy & ratio$CopyNumber!= -1)
      points(ratio$Start[tt],ratio$Ratio[tt]*ploidy,pch = ".",col = colors()[461])
      tt <- which(ratio$Chromosome==i)
      
      #UNCOMMENT HERE TO SEE THE PREDICTED COPY NUMBER LEVEL:
      #points(ratio$Start[tt],ratio$CopyNumber[tt], pch = ".", col = colors()[24],cex=4)
      
    }
    tt <- which(ratio$Chromosome==i)
    
    #UNCOMMENT HERE TO SEE THE EVALUATED MEDIAN LEVEL PER SEGMENT:
    #points(ratio$Start[tt],ratio$MedianRatio[tt]*ploidy, pch = ".", col = colors()[463],cex=4)
    
  }  
  dev.off()  
  
} else {cat ("WARNING: To get a .png image with copy number profile, you can provide as input a file with suffix 'ratio.txt'\n")}


#------------------------------------------------------

#plot .png for the _BAF.txt file:


if (BAFfileInd) {
  BAF <-read.table(args[BAFfileInd], header=TRUE);
	BAF<-data.frame(BAF)

	png(filename = paste(args[BAFfileInd],".png",sep = ""), width = 1180, height = 1180,
	    units = "px", pointsize = 20, bg = "white", res = NA)
	plot(1:10)
	op <- par(mfrow = c(5,5))

	for (i in c(1:22,'X','Y')) {
	    tt <- which(BAF$Chromosome==i)
	    if (length(tt)>0){
		lBAF <-BAF[tt,]
		plot(lBAF$Position,lBAF$BAF,ylim = c(-0.1,1.1),xlab = paste ("position, chr",i),ylab = "BAF",pch = ".",col = colors()[1])

		tt <- which(lBAF$A==0.5)		
		points(lBAF$Position[tt],lBAF$BAF[tt],pch = ".",col = colors()[92])
		tt <- which(lBAF$A!=0.5 & lBAF$A>=0)
		points(lBAF$Position[tt],lBAF$BAF[tt],pch = ".",col = colors()[62])
		tt <- 1
		pres <- 1

		if (length(lBAF$A)>4) {
			for (j in c(2:(length(lBAF$A)-pres-1))) {
				if (lBAF$A[j]==lBAF$A[j+pres]) {	
					tt[length(tt)+1] <- j 
				}
			}
			points(lBAF$Position[tt],lBAF$A[tt],pch = ".",col = colors()[24],cex=4)
			points(lBAF$Position[tt],lBAF$B[tt],pch = ".",col = colors()[24],cex=4)	
		}

		tt <- 1
		pres <- 1
		if (length(lBAF$FittedA)>4) {
			for (j in c(2:(length(lBAF$FittedA)-pres-1))) {
				if (lBAF$FittedA[j]==lBAF$FittedA[j+pres]) {	
					tt[length(tt)+1] <- j 
				}
			}
			points(lBAF$Position[tt],lBAF$FittedA[tt],pch = ".",col = colors()[463],cex=4)
			points(lBAF$Position[tt],lBAF$FittedB[tt],pch = ".",col = colors()[463],cex=4)	
		}

	   }

	}
	dev.off()

} else {cat ("WARNING: To get a .png image with BAF profile, you can provide as input a file with suffix 'BAF.txt'\n")}
