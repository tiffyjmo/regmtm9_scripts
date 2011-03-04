#written by Tiffany Morris
#June 2010

setwd(wd)

chromosome = c('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22')

for (chr in chromosome) 
{ 
  ampsCounter = 0
  delCounter = 0
  gainCounter = 0
  lossCounter = 0
  
  workingDirectory<-paste(wd,"/Regions/",sep="")
  setwd(workingDirectory) 

  theFilters = c("chromosome_name", "start", "end")
  theAttributes = c("mirbase_id")
    
  ampRegionData<-paste("computedAmpRegions_",chr,".txt",sep="")
  ampRegionData<-read.table(ampRegionData,header=TRUE,stringsAsFactors=FALSE,sep="\t")

  delRegionData<-paste("computedDelRegions_",chr,".txt",sep="")
  delRegionData<-read.table(delRegionData,header=TRUE,stringsAsFactors=FALSE,sep="\t")

  gainRegionData<-paste("computedGainRegions_",chr,".txt",sep="")
  gainRegionData<-read.table(gainRegionData,header=TRUE,stringsAsFactors=FALSE,sep="\t")

  lossRegionData<-paste("computedLossRegions_",chr,".txt",sep="")
  lossRegionData<-read.table(lossRegionData,header=TRUE,stringsAsFactors=FALSE,sep="\t")

  workingDirectory<-paste(wd,"/mirnaLists/",sep="")
  library(biomaRt)

  ensembl = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  
  if(nrow(ampRegionData) != 0) 
  {
    ampGenes <- getBM(attributes = theAttributes, filters = theFilters, values = list(ampRegionData$CHR,ampRegionData$STARTPOS,ampRegionData$ENDPOS), mart=ensembl)
    #ampGenes=ampGenes[order(ampGenes$chromosome_name,ampGenes$start_position),]
    file1 <- paste(workingDirectory,study,"_AmpsMiRNAs_",chr,".txt",sep="")
    sort.ampGenes<-as.data.frame(ampGenes[order(ampGenes$mirbase_id) , ])
    colnames(sort.ampGenes)<-"mirbase_id"
    write.table(sort.ampGenes, file=file1, row.name=FALSE, quote=FALSE, sep="\t")
    if(ampsCounter==0)
    {
      mergedAmps = ampGenes
      ampsCounter = ampsCounter +1
    }
    else
    {
      mergedAmps = merge(mergedAmps,ampGenes,by=c("mirbase_id"),all=TRUE)
    }
  }
  if(nrow(gainRegionData) != 0)
  {
    gainGenes <- getBM(attributes = theAttributes, filters = theFilters, values = list(gainRegionData$CHR,gainRegionData$STARTPOS,gainRegionData$ENDPOS), mart=ensembl)
    #gainGenes=gainGenes[order(gainGenes$chromosome_name,gainGenes$start_position),]
    file2 <- paste(workingDirectory,study,"_GainMiRNAs_",chr,".txt",sep="")
    sort.gainGenes<-as.data.frame(gainGenes[order(gainGenes$mirbase_id) , ])
    colnames(sort.gainGenes)<-"mirbase_id"
    write.table(sort.gainGenes, file=file2, row.name=FALSE, quote=FALSE, sep="\t")
    if(gainCounter==0)
    {
      mergedGain = gainGenes
      gainCounter = gainCounter +1
    }
    else
    {
      mergedGain = merge(mergedGain,gainGenes,by=c("mirbase_id"),all=TRUE)
    }
  }
  if(nrow(delRegionData) != 0)
  {
    delGenes <- getBM(attributes = theAttributes, filters = theFilters, values = list(delRegionData$CHR,delRegionData$STARTPOS,delRegionData$ENDPOS), mart=ensembl)
    #delGenes=delGenes[order(delGenes$chromosome_name,delGenes$start_position),]
    file3 <- paste(workingDirectory,study,"_DelMiRNAs_",chr,".txt",sep="")
    sort.delGenes<-as.data.frame(delGenes[order(delGenes$mirbase_id) , ])
    colnames(sort.delGenes)<-"mirbase_id"
    write.table(sort.delGenes, file=file3, row.name=FALSE, quote=FALSE, sep="\t")
    if(delCounter==0)
    {
      mergedDel = delGenes
      delCounter = delCounter +1
    }
    else
    {
      mergedDel = merge(mergedDel,delGenes,by=c("mirbase_id"),all=TRUE)
    }
  }
  if(nrow(lossRegionData) != 0)
  {
    lossGenes <- getBM(attributes = theAttributes, filters = theFilters, values = list(lossRegionData$CHR,lossRegionData$STARTPOS,lossRegionData$ENDPOS), mart=ensembl)
    #lossGenes=lossGenes[order(lossGenes$chromosome_name,lossGenes$start_position),]
    file4 <- paste(workingDirectory,study,"_LossMiRNAs_",chr,".txt",sep="")
    sort.lossGenes<-as.data.frame(lossGenes[order(lossGenes$mirbase_id) , ])
    colnames(sort.lossGenes)<-"mirbase_id"
    write.table(sort.lossGenes, file=file4, row.name=FALSE, quote=FALSE, sep="\t")
    if(lossCounter==0)
    {
      mergedLoss = lossGenes
      lossCounter = lossCounter +1
    }
    else
    {
      mergedLoss = merge(mergedLoss,lossGenes,by=c("mirbase_id"),all=TRUE)
    }
  }
  print(study)
  print(chr)
  print("Done")
  
  file5 <- paste(workingDirectory,"mergedMiRNAAmps_",chr,".txt",sep="")
#  mergedAmps<-mergedAmps[order(mergedAmps$hgnc_symbol) , ]
  write.table(mergedAmps,file=file5, row.name=FALSE, quote=FALSE, sep="\t")
  file6 <- paste(workingDirectory,"mergedMiRNAGain_",chr,".txt",sep="")
#  mergedGain<-mergedGain[order(mergedGain$hgnc_symbol) , ]
  write.table(mergedGain,file=file6, row.name=FALSE, quote=FALSE, sep="\t")
  file7 <- paste(workingDirectory,"mergedMiRNADel_",chr,".txt",sep="")
#  mergedDel<-mergedDel[order(mergedDel$hgnc_symbol) , ]
  write.table(mergedDel,file=file7, row.name=FALSE, quote=FALSE, sep="\t")
  file8 <- paste(workingDirectory,"mergedMiRNALoss_",chr,".txt",sep="")
#  mergedLoss<-mergedLoss[order(mergedLoss$hgnc_symbol) , ]
  write.table(mergedLoss,file=file8, row.name=FALSE, quote=FALSE, sep="\t")
}
