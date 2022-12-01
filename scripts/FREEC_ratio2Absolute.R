args <- commandArgs()

#read _ratio.txt:
dataTable <-read.table(args[4], header=TRUE);

ratio<-data.frame(dataTable)

cat ("Chromosome\tStart\tEnd\tNum_Probes\tSegment_Mean")

weirdNum=345

for (i in c(unique(dataTable$Chromosome))) {
	tt <- which(ratio$Chromosome==i)
	if (length(tt)>0) {
	  myMedianRatio=weirdNum
	  mySegNUmber=0
	  for (j in c(1:length(tt))) {
	    if (ratio$MedianRatio[tt][j]!=myMedianRatio) {
	      #print the record if it is not -1 or weirdNum
	      if (myMedianRatio != -1 & myMedianRatio != weirdNum) {
	        cat(i,"\t", ratio$Start[tt][j-mySegNumber], "\t", ratio$Start[tt][j]-1,"\t",mySegNumber,"\t",log(myMedianRatio,2),"\n")
	      }
	      #update
	      mySegNumber=1
	      myMedianRatio=ratio$MedianRatio[tt][j]
	    }
	    else {
	      #update:
	      mySegNumber = mySegNumber+1
	    }
	  }
	  #print the last segment
	  if (myMedianRatio != -1 & myMedianRatio != weirdNum) {
	    cat(i,"\t", ratio$Start[tt][j-mySegNumber], "\t", ratio$Start[tt][j-1],"\t",mySegNumber,"\t",log(myMedianRatio,2),"\n")
	  }
	}	  
}
