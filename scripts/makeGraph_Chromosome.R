args <- commandArgs()

dataTable <-read.table(args[6], header=TRUE);
print( paste (args[6],"read"))
ratio<-data.frame(dataTable)
chr <- type.convert(args[4])
#chr <- 'X'
ploidy <- type.convert(args[5])

maxLevelToPlot <- 3
for (i in c(1:length(ratio$Ratio))) {
	if (ratio$Ratio[i]>maxLevelToPlot) {
		ratio$Ratio[i]=maxLevelToPlot;
	}
}

png(filename = paste(args[6],".png",sep = ""), width = 1180, height = 1180,
    units = "px", pointsize = 20, bg = "white", res = NA)
plot(1:10)
op <- par(mfrow = c(2,1))
i <- chr
	tt <- which(ratio$Chromosome == i)
	if (length(tt)>0) {
		 plot(ratio$Start[tt],ratio$Ratio[tt]*ploidy,ylim = c(0,maxLevelToPlot*ploidy),xlab = paste ("position, chr",i),ylab = "normalized copy number profile",pch = ".",col = colors()[88],cex=2)
		 tt <- which(ratio$Chromosome==i  & ratio$CopyNumber>ploidy )
		 points(ratio$Start[tt],ratio$Ratio[tt]*ploidy,pch = ".",col = colors()[136],cex=2)
		 tt <- which(ratio$Chromosome==i  & ratio$Ratio==maxLevelToPlot)
		 points(ratio$Start[tt],ratio$Ratio[tt]*ploidy,pch = ".",col = colors()[136],cex=4)
		 tt <- which(ratio$Chromosome==i  & ratio$CopyNumber<ploidy )
		 points(ratio$Start[tt],ratio$Ratio[tt]*ploidy,pch = ".",col = colors()[461],cex=2)
		 tt <- which(ratio$Chromosome==i)
		 points(ratio$Start[tt],ratio$CopyNumber[tt], pch = ".", col = colors()[24],cex=2)
		 points(ratio$Start[tt],ratio$MedianRatio[tt]*ploidy, pch = ".", col = colors()[98],cex=4)
	}
if (length(args)>=7) {
	dataTable <-read.table(args[7], header=TRUE);
	BAF<-data.frame(dataTable)
	tt <- which(BAF$Chromosome==i)
	lBAF <-BAF[tt,]
	plot(lBAF$Position,lBAF$BAF,ylim = c(-0.1,1.1),xlab = paste ("position, chr",i),ylab = "BAF",pch = ".",col = colors()[1])
	tt <- which(lBAF$A==0.5)		
	points(lBAF$Position[tt],lBAF$BAF[tt],pch = ".",col = colors()[92])
	tt <- which(lBAF$A!=0.5 & lBAF$A>=0)
	points(lBAF$Position[tt],lBAF$BAF[tt],pch = ".",col = colors()[450])
	tt <- which(lBAF$B==1)
	points(lBAF$Position[tt],lBAF$BAF[tt],pch = ".",col = colors()[62])
	tt <- 1
	pres <- 1
	for (j in c(2:(length(lBAF$A)-pres-1))) {
			if (lBAF$A[j]==lBAF$A[j+pres]) {	
				tt[length(tt)+1] <- j 
			}
	}
	points(lBAF$Position[tt],lBAF$A[tt],pch = ".",col = colors()[24],cex=4)
	points(lBAF$Position[tt],lBAF$B[tt],pch = ".",col = colors()[24],cex=4)
	print( paste (args[7],"read"))
}
dev.off()

