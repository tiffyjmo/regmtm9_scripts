#written by Claude Chelala
#adapted by Tiffany Morris April 2010

setwd(wd)

#par(mfrow=c(4,1));

patientFile<-paste(wd,"/freq_patient.txt",sep="")
cellFile<-paste(wd,"/freq_cell.txt",sep="")


patientData <- read.table(patientFile, sep="\t", header = TRUE, as.is = TRUE)
patientSubset = subset(patientData, select = c('probeID','chr','Pos','patient_Freq_Gain','patient_Freq_Loss'))
patientValues <- patientData[,c(1:(numPatients+3))]

cellData <- read.table(cellFile, sep="\t", header = TRUE, as.is = TRUE)
cellSubset = subset(cellData, select = c("probeID","chr","Pos","cell_Freq_Gain","cell_Freq_Loss"))
cellValues <- cellData[,c(1:(numCell+3))]

freqData = merge(patientSubset,cellSubset,by=c("probeID","chr","Pos"),all=TRUE)

#all data (with log2 ratios)                                                             
Allvalues = merge(patientValues,cellValues,by=c("probeID","chr","Pos"),all=TRUE)
All = merge(Allvalues,freqData,by=c("probeID","chr","Pos"),all=TRUE)
All = All[order(All$chr,All$Pos),]

values = Allvalues[,c(4:(numPatients+numCell+3))]

file <- paste(wd,"/AllFreq.txt",sep="")
write.table(All, file=file, row.name=FALSE, quote=FALSE, sep="\t")

file2 <- paste(wd,"/",study,"_log2.txt",sep="")
write.table(Allvalues, file=file2, row.name=FALSE, quote=FALSE, sep="\t")

###############FREQ PLOT##########################
#label chromosomes
png("FreqPlot.jpeg")
#pdf("FreqPlot.pdf"

#needs regions file!!
freqFile<-paste(wd,"/Regions/filteredRegionsFreqPlot.txt",sep="")

data = All
data <- read.table(freqFile, sep="\t", header = TRUE, as.is = TRUE)
data = data[order(data$chr,data$Pos),]
labels_chr <- data.matrix(summary(as.factor(data$chr)))

test1<- data.frame(labels_chr,row.names(labels_chr) )
test <- data.frame(unique(data$chr))
colnames(test) = c("chr")
colnames(test1) = c("count","chr")
F1 <- merge(test,test1, by="chr", sort=FALSE)
for(i in 2:length(row.names(F1))){F1[i,2] = F1[i-1,2] + F1[i,2] ; }

F1$label <- NULL ; F1[1,3] <- F1[1,2] / 2 ;
for (i in 2:length(row.names(F1))){ F1[i,3] <- (F1[i,2]+F1[i-1,2])/2; }

y1=100*((data$patient_Freq_Gain+data$cell_Freq_Gain)/samples)
y2=paste("-",100*((data$patient_Freq_Loss+data$cell_Freq_Loss)/samples),sep="")

setwd(wd)
#pdf("FreqPlot.pdf", horizontal =TRUE,width = 10.0, height = 9.0)

#plot gain
graphtitle <- paste("Frequency Plot of Genome Study ",platform," ",study,sep="")
plot(y1, type='h',  xaxt="n",  yaxt="n", col="green", main = graphtitle , ylim=range(-100, 100), xlab='Chromosome Number',  ylab='Fraction of Patient for Gain or Loss',xaxs = "i", yaxs = "i")

#plot loss
points(y2, type='h', col="red")

#label for chromosomes
x= c(-100,-80,-60,-40,-20,0,20,40,60,80,100)
y= c(1:length(F1[,2]))
axis(1, at = c(F1[,2]), label =FALSE, tick = TRUE, las = 1, col = "black", lty = "dotted", tck = 1 );
axis(1, at = c(F1[,3]), label =F1$chr, tick = FALSE );
axis(2, at = c(x), labels=x, tick = TRUE, las = 1, col = "black", lty = "dotted", tck = 1 );

dev.off()  

###############AMP/DEL FREQ PLOT##########################
#label chromosomes
png("AMPDEL_FreqPlot.jpeg")
#pdf("AMPDEL_FreqPlot.pdf")
labels_chr <- data.matrix(summary(as.factor(ampDelValues$chr)))

test1<- data.frame(labels_chr,row.names(labels_chr) )
test <- data.frame(unique(ampDelValues$chr))
colnames(test) = c("chr")
colnames(test1) = c("count","chr")
F1 <- merge(test,test1, by="chr", sort=FALSE)
for(i in 2:length(row.names(F1))){F1[i,2] = F1[i-1,2] + F1[i,2] ; }

F1$label <- NULL ; F1[1,3] <- F1[1,2] / 2 ;
for (i in 2:length(row.names(F1))){ F1[i,3] <- (F1[i,2]+F1[i-1,2])/2; }

y1=100*((ampDelValues$patient_Freq_Amp+ampDelValues$cell_Freq_Amp)/samples)
y2=paste("-",100*((ampDelValues$patient_Freq_Del+ampDelValues$cell_Freq_Del)/samples),sep="")

setwd(wd)
#pdf("FreqPlot.pdf", horizontal =TRUE,width = 10.0, height = 9.0)

#plot gain
graphtitle <- paste("Amp/Del Frequency Plot of Genome Study ",study,sep="")
plot(y1, type='h',  xaxt="n",  yaxt="n", col="green", main = graphtitle , ylim=range(-100, 100), xlab='Chromosome Number',  ylab='Fraction of Patient for Gain or Loss',xaxs = "i", yaxs = "i")

#plot loss
points(y2, type='h', col="red")

#label for chromosomes
x= c(-100,-80,-60,-40,-20,0,20,40,60,80,100)
y= c(1:length(F1[,2]))
abline(h=0.2,col="grey")
axis(1, at = c(F1[,2]), label =FALSE, tick = TRUE, las = 1, col = "black", lty = "dotted", tck = 1 );
axis(1, at = c(F1[,3]), label =F1$chr, tick = FALSE );
axis(2, at = c(x), labels=x, tick = TRUE, las = 1, col = "black", lty = "dotted", tck = 1 );

dev.off() 


  
###############FREQ PLOT Per Chromosome##########################
setwd(wd)
png("perChrom.jpeg");
#pdf("perChrom.pdf");
par(mfrow=c(4,1));
#par(mfrow=c(2,2));
#chromosomes =c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22)
chromosomes =c(F1$chr)
for(test in chromosomes) 
{
  data.plot.subset = subset(data, chr == test)
  data.plot.subset = data.plot.subset[order(data.plot.subset$chr,data.plot.subset$Pos),]

  y1=100*((data.plot.subset$patient_Freq_Gain+data.plot.subset$cell_Freq_Gain)/samples)
  y2=paste("-",100*((data.plot.subset$patient_Freq_Loss+data.plot.subset$cell_Freq_Loss)/samples),sep="")
  #postscript("Freq_PET_all_chr.ps", horizontal =TRUE,width = 10.0, height = 9.0)

  #plot gain
  graphtitle <- paste("Frequency Plot per Chromosome of Genome Study ",study,sep="")
  plot(y1, type='h',  xaxt="n",  yaxt="n", col="green", main = graphtitle, xlim=range(0, 200),ylim=range(-100, 100), xlab=test,  ylab='Fraction of Patient for Gain or Loss',xaxs = "i", yaxs = "i")

  #plot loss
  points(y2, type='h', col="red")

  #label for chromosomes
  x= c(-100,-80,-60,-40,-20,0,20,40,60,80,100)
  y = c(1)
  axis(2, at = c(x), labels=x, tick = TRUE, las = 1, col = "black", lty = "dotted", tck = 1 );
}

dev.off()
#######################Clustering with A2R#################################################
