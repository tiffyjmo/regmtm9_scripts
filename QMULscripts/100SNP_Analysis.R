#Tiffany Morris August 2010

##############100KSNP analysis###############################################
setwd(analysisDir)
library("aroma.affymetrix")
log <- Arguments$getVerbose(-8, timestamp=TRUE)
options(digits=4)

#check platform type
chipTypes <- c("Mapping50K_Hind240","Mapping50K_Xba240")

#####check files are all in place
#####change working directory
cdfs <- lapply(chipTypes, FUN=function(chipType) {
  AffymetrixCdfFile$byChipType(chipType)
})
print(cdfs)

gis <- lapply(cdfs, getGenomeInformation)
print(gis)

sis <- lapply(cdfs, getSnpInformation)
print(sis)

#####Defining CEL set
cesNList <-  list()

chipType <- chipTypes[1] #repeat with 2

cs <- AffymetrixCelSet$byName(dataSetName, chipType=chipType)

csR <- extract(cs, !isDuplicated(cs))

print(cs)

#####Quality assessment
#cs<-csR
#par(mar=c(4,4,1,1)+0.1)
#plotDensity(cs, lwd=2, ylim=c(0,0.40))
#stext(side=3, pos=0, getFullName(cs))

######Allelic crosstalk calibration
acc <- AllelicCrosstalkCalibration(csR, model="CRMAv2")
print(acc)

csC <- process(acc, verbose=verbose)
print(csC)

######Quality assessment post-normalisation
#cs <- csC
#par(mar=c(4,4,1,1)+0.1)
#plotDensity(cs, lwd=2, ylim=c(0,0.40))
#stext(side=3, pos=0, getFullName(cs))

#####Normalisation for nucleotide-position probe sequence effects
bpn <- BasePositionNormalization(csC, target="zero")  
print(bpn)

csN <- process(bpn, verbose=verbose) #takes 2-8min/array
print(csN)

#####Quality Assessment	post-normalisation2
#cs <- csN
#par(mar=c(4,4,1,1)+0.1)
#plotDensity(cs, lwd=2, ylim=c(0,0.40))
#stext(side=3, pos=0, getFullName(cs))


#####Probe Summarisation
plm <- RmaCnPlm(csN, mergeStrands=TRUE, combineAlleles=TRUE) #RmaCnPlm is more appropriate for 100k chip
print(plm)
 
fit(plm, verbose=verbose) #below code should be faster??

#if (length(findUnitsTodo(plm)) > 0) {
  # Fit CN probes quickly (~5-10s/array + some overhead)
#  units <- fitCnProbes(plm, verbose=verbose)
#  str(units)
  # int [1:945826] 935590 935591 935592 935593 935594 935595 ...
  
  # Fit remaining units, i.e. SNPs (~5-10min/array)
#  units <- fit(plm, verbose=verbose)
#  str(units)
#}
  
ces <- getChipEffectSet(plm)
print(ces)

#####Normalisation for PCR fragment-length effects
fln <- FragmentLengthNormalization(ces, target="zero") 
print(fln)


cesNList[[chipType]] <- process(fln, verbose=verbose)

#use this if combineAlleles is set to FALSE
#asN <- AlleleSummation(cesNList[[chipType]])
#cesNTotal[[chipType]] <- process(asN, verbose=log)

######################################################Repeat above for second chip######################################
chipType <- chipTypes[2] #repeat with 2

cs <- AffymetrixCelSet$byName(dataSetName, chipType=chipType)
csR <- extract(cs, !isDuplicated(cs))
print(cs)

#####Quality assessment
#cs<-csR
#par(mar=c(4,4,1,1)+0.1)
#plotDensity(cs, lwd=2, ylim=c(0,0.40))
#stext(side=3, pos=0, getFullName(cs))

######Allelic crosstalk calibration
acc <- AllelicCrosstalkCalibration(csR, model="CRMAv2")
print(acc)

csC <- process(acc, verbose=verbose)
print(csC)

######Quality assessment post-normalisation
#cs <- csC
#par(mar=c(4,4,1,1)+0.1)
#plotDensity(cs, lwd=2, ylim=c(0,0.40))
#stext(side=3, pos=0, getFullName(cs))

#####Normalisation for nucleotide-position probe sequence effects
bpn <- BasePositionNormalization(csC, target="zero")
print(bpn)

csN <- process(bpn, verbose=verbose) #takes 2-8min/array
print(csN)

#####Quality Assessment	
#cs <- csN
#par(mar=c(4,4,1,1)+0.1)
#plotDensity(cs, lwd=2, ylim=c(0,0.40))
#stext(side=3, pos=0, getFullName(cs))


#####Probe Summarisation
plm <- RmaCnPlm(csN, mergeStrands=TRUE, combineAlleles=TRUE) #RmaCnPlm is more appropriate for 100k chip
print(plm)
 
fit(plm, verbose=verbose) #below code should be faster??

#if (length(findUnitsTodo(plm)) > 0) {
  # Fit CN probes quickly (~5-10s/array + some overhead)
#  units <- fitCnProbes(plm, verbose=verbose)
#  str(units)
  # int [1:945826] 935590 935591 935592 935593 935594 935595 ...
  
  # Fit remaining units, i.e. SNPs (~5-10min/array)
#  units <- fit(plm, verbose=verbose)
#  str(units)
#}
  
ces <- getChipEffectSet(plm)
print(ces)

#####Normalisation for PCR fragment-length effects
fln <- FragmentLengthNormalization(ces)  #target="zero" not sure appropriate for this platform???
print(fln)


cesNList[[chipType]] <- process(fln, verbose=verbose)

#use this if combineAlleles is set to FALSE
#asN <- AlleleSummation(cesNList[[chipType]])
#cesNTotal[[chipType]] <- process(asN, verbose=log)

##########################################################END REPEAT##################################
print(cesNList)

#use this if combineAlleles is set to FALSE
#cesNList<-cesNTotal
 
#####Calculation of Raw Copy Numbers using a common reference

idxnorm <- (1:totalNormals);
idxtumours <- ((totalNormals+1):(totalNormals+totalSamples));

cesN1<-cesNList[[1]]
norms1<-extract(cesN1, idxnorm);

cesN2<-cesNList[[2]]
norms2<-extract(cesN2, idxnorm);

tumours1<-extract(cesN1, idxtumours);
tumours2<-extract(cesN2, idxtumours);
 
########Fit model
#avgnorm <- getAverageFile(norms1, verbose=log);
#cbs <- CbsModel(tumours1,avgnorm);
#rawCNs <- extractRawCopyNumbers(cbs, array=1, chromosome=1)
#rawCNs <- as.data.frame(rawCNs)
#str(rawCNs)


#######for 1st array...extract raw copy numbers
#check platform type again??
cdf1 <- AffymetrixCdfFile$byChipType("Mapping50K_Hind240")
cdf2 <- AffymetrixCdfFile$byChipType("Mapping50K_Xba240")

gi1 <- getGenomeInformation(cdf1)
gi2 <- getGenomeInformation(cdf2)

si1 <- getSnpInformation(cdf1)
si2 <- getSnpInformation(cdf2)

#####loop through arrays to compute Log2Ratios
for(array in 1:totalSamples)
{
  allt={}
  currentTumour1 <- getFile(tumours1,array)
  name1 = getName(currentTumour1)
  currentTumour2 <- getFile(tumours2,array)
  name2 = getName(currentTumour2)
  name <- paste(name1,"_",name2,sep="")

  for(Chromosome in 1:22) 
  {

	   units <- getUnitsOnChromosome(gi1, chromosome=Chromosome)
	   Position <- getPositions(gi1, units=units)
	   ProbeID <- getUnitNames(cdf1, units=units)
	   theta <- extractMatrix(currentTumour1, units=units)
     thetaR <- extractMatrix(norms1, units=units)
	   thetaR <- rowMedians(thetaR, na.rm=TRUE)

	   M = log2(theta/thetaR) 

	   LogRatio <- M[,1]
	   array1 <- cbind(ProbeID,Chromosome,Position,LogRatio)
	   array1<-as.data.frame(array1)
     array1$Chromosome<-as.character(array1$Chromosome)	
     array1$Chromosome<-as.numeric(array1$Chromosome)
     array1$Position<-as.character(array1$Position)	
     array1$Position<-as.numeric(array1$Position) 

	   ####repeat for 2nd array

	   units <- getUnitsOnChromosome(gi2, chromosome=Chromosome)
	   Position <- getPositions(gi2, units=units)
	   ProbeID <- getUnitNames(cdf2, units=units)
	   theta <- extractMatrix(currentTumour2, units=units)
           thetaR <- extractMatrix(norms2, units=units)
	   thetaR <- rowMedians(thetaR, na.rm=TRUE)

	   M = log2(theta/thetaR) 

	   LogRatio <- M[,1]
	   array2 <- cbind(ProbeID,Chromosome,Position,LogRatio)
	   array2<-as.data.frame(array2)
     array2$Chromosome<-as.character(array2$Chromosome)	
     array2$Chromosome<-as.numeric(array2$Chromosome)
     array2$Position<-as.character(array2$Position)	
     array2$Position<-as.numeric(array2$Position) 
	   
	
	   #####combine data for two arrays
	   allt <- rbind(allt,array1,array2)
    
  }
  #allt<-allt[order(as.numeric(allt$Chromosome),as.numeric(allt$Position)),]
	file <- paste(resultsDir,"/CGHweb_InputFiles/",name,"_Input",".txt",sep="")
  write.table(allt,file=file,append=TRUE,sep="\t",row.name=FALSE,quote=FALSE)
}
