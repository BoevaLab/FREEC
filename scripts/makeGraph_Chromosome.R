path = "/data/tmp/cgurjao/Graphical_outputs/R\ scripts/FREEC_7.4/Single-end"
setwd(path)

args = c()
args[8] = 0
args[7] = 0
args[6] = "ratios_0.txt"
args[5] = 2
args[4] = 1

dataTable <-read.table(args[6], header=T, na.strings="NA");
print(paste (args[6],"read"))
ratio<-data.frame(dataTable)
chr <- type.convert(args[4])
ploidy <- type.convert(args[5])
CN_subc = ratio$Subclones

max_value = 10#max(ratio$CopyNumber)

plot(1:10)
even_odd = 0
even_odd2 = 0
par(mfrow = c(2,1))
for (i in c(1:11))
{
  png(filename = paste(path, "/Results/","chr", i,".png",sep = ""), width = 1180, height = 1180,
    units = "px", pointsize = 20, bg = "white", res = NA)
	tt <- which(ratio$Chromosome==i)
	if (length(tt)>0) {
	   tt <- which(ratio$Chromosome==i  & ratio$CopyNumber==ploidy )
		 plot(ratio$Start[tt],ratio$Ratio[tt]*ploidy,ylim = c(0,max_value),xlab = paste ("position, chr",i),ylab = "normalized copy number profile",pch = ".",col = rgb(0,1,0, alpha = 0.2),cex=8)
		 tt <- which(ratio$Chromosome==i  & ratio$CopyNumber>ploidy )
		 points(ratio$Start[tt],ratio$Ratio[tt]*ploidy,pch = ".",col = rgb(1,0,0, alpha = 0.2),cex=8)
		 tt <- which(ratio$Chromosome==i  & ratio$CopyNumber<ploidy )
		 points(ratio$Start[tt],ratio$Ratio[tt]*ploidy,pch = ".",col = rgb(0,0,1, alpha = 0.2),cex=8)
		 tt <- which(ratio$Chromosome==i)
		 points(ratio$Start[tt],ratio$CopyNumber[tt], pch = ".", col = colors()[24],cex=8)
		 tt <- which(ratio$Chromosome==i)
		 for (k in c(1:length(levels(CN_subc))))
		 {
		   if (levels(CN_subc)[k] != "0/0" && length(CN_subc) > 0)
		   {
		     ttt <- which(ratio$Subclones == levels(CN_subc)[k] &  ratio$Chromosome==i)
         if (length(ttt)>0)
		        {
		     j = 1
		     CN = ''
         pop = ''
		     while (substr(as.character(ratio$Subclones[ttt]),j,j)[1] != "/")
		     { 
		       CN = paste(CN,substr(as.character(ratio$Subclones[ttt]),j,j))
		       j = j+1
		     }
         j = 9
		     while (substr(as.character(ratio$Subclones[ttt]),j,j)[1] != "/" )
		     { 
		        pop[j] = paste(pop,substr(as.character(ratio$Subclones[ttt]),j,j))
            j = j-1
		     }

         x = NA
         y = NA
         if (round(length(ttt)/2) > 0 && (-ratio$Start[ttt[1]] + ratio$Start[ttt[length(ttt)-1]]) > 10000000 )
           {
           if (even_odd %% 2 == 0) #(ratio$CopyNumber[ttt][round(length(ttt)/2)] > as.numeric(gsub(" ", "",CN[1])))
		        {
             points(ratio$Start[ttt],CN, pch = ".", col = colors()[99],cex=8)
             if (even_odd2 %% 2 == 0)
                {x = as.numeric(ratio$Start[ttt][round(length(ttt)/3)])
                y = as.numeric(CN[1]) - 0.25
                even_odd2 = even_odd2 + 1}
             else
                {x = as.numeric(ratio$Start[ttt][round(length(ttt)/3)])
                 y = as.numeric(CN[1]) + 0.25
                 even_odd2 = even_odd2 + 1}
            even_odd = even_odd +1   
            }
           else 
            {
              points(ratio$Start[ttt],CN, pch = ".", col = colors()[96],cex=8)
              if (even_odd2 %% 2 != 0)
                {
                x = as.numeric(ratio$Start[ttt][round(length(ttt)/3)])
                y = as.numeric(CN[1]) + 0.25
                even_odd2 = even_odd2 + 1}
              else
                {x = as.numeric(ratio$Start[ttt][round(length(ttt)/3)])
                y = as.numeric(CN[1]) - 0.25
                even_odd2 = even_odd2 + 1}
            even_odd = even_odd +1  
            } 
          }
         pop = as.numeric(gsub("NA", "",gsub(",","",gsub(' ',"", toString(pop))))) * 100;
         pop = paste(as.character(as.integer(pop)),"%")
         text(x, y, as.character(pop), cex = 0.8, col = 'black')
        }
      }
		 }
   }             
             
  if (args[7] != 0) {
  dataTable <-read.table(args[7], header=TRUE);
	BAF<-data.frame(dataTable)
	tt <- which(BAF$Chromosome==i)
	lBAF <-BAF[tt,]
	plot(lBAF$Position,lBAF$BAF, ylim = c(-0.1,1.1),xlab = paste ("position, chr",i),ylab = "BAF",pch = ".",col = colors()[1])
	tt <- which(lBAF$A==0.5)		
	points(lBAF$Position[tt],lBAF$BAF[tt],pch = ".",col = colors()[92])
	tt <- which(lBAF$A!=0.5 & lBAF$A>=0)
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
}
