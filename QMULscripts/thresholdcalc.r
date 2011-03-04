#written by Claude

setwd(wd); #wd is the specific results directory for the current study

#--data input--
file<-paste(wd,"/allFreq.txt",sep="")

print("calculating threshold...")

#--read table--
data <- read.table(file, sep="\t", header = TRUE)

k = stack(data[,first:columns])#first to last patient
rowmed = k$values
q3 =  quantile ( rowmed, .75, na.rm = T) #third quantile
q1 = quantile ( rowmed, .25, na.rm = T)   # 1st quantile
  
iqr = rowmed[rowmed>q1 & rowmed < q3]
length(iqr)
gain_loss =  3* sd(iqr,na.rm = TRUE)

High = rowmed[rowmed > gain_loss  | rowmed < -gain_loss];
amp =  quantile ( High, .97, na.rm = T); #third quantile
del = quantile ( High, .03, na.rm = T);   # 1st quartile
quantile(High, probs = seq(0, 1, 0.01), na.rm = T, names = TRUE);

setwd(wd)
png("density.jpeg");
#pdf("density.pdf");
plot(density(rowmed, na.rm = T,adjust=0.5));

abline( v =  gain_loss, col = 'green', lwd = 2);
abline( v = -gain_loss, col = 'red', lwd = 2);

#legend( locator(1), legend = c("gain","loss"),
#        lty = c(1,1) ,    
#        lwd = 2, col = c('blue','red'))

abline( v =  gain_loss, col = 'green', lwd = 2)
abline( v = -gain_loss, col = 'red', lwd = 2)

#legend( locator(1), legend = c("gain","loss"),
#	lty = c(1,1) ,    
#        lwd = 2, col = c('blue','red'))

abline( v = quantile( High, probs = 0.03, na.rm = T), col = 'red', lwd = 2)
abline( v = quantile( High, probs = 0.97, na.rm = T), col = 'green', lwd = 2)
  
#legend( locator(1), legend = c("Amplification","Deletion"),
#         lty = c(1,1) ,    # line type, 1 is regular
#         lwd = 2, col = c('blue','red'))

dev.off()

pdf("histogram.pdf")
hist(rowmed, breaks=30) #
dev.off()
pdf("histogram2.pdf")
hist(rowmed, freq=F, density=10, angle=45, col="blue")
abline( v = gain_loss, col = 'green', lwd = 2)
abline( v = -gain_loss, col = 'red', lwd = 2)
abline( v = quantile( High, probs = 0.03,na.rm=T), col = 'red', lwd = 2)
abline( v = quantile( High, probs = 0.97,na.rm=T), col = 'green', lwd = 2)
dev.off()

#send this information to the report
IQR<-length(iqr) ;sd_IQR<-sd(iqr,na.rm = TRUE);
message2<-paste("The length of the IQR is: ",IQR,"The SD of the IQR is: ",sd_IQR)
print(message2)


print (message)
print ("Quantile_rowmed")
print(quantile(rowmed, probs = seq(0, 1, 0.01), na.rm = TRUE, names = TRUE))

gain_loss = 3*sd(iqr,na.rm = TRUE)
message<-paste("The threshold is: ",gain_loss)
min = min(rowmed, na.rm=T)
max = max(rowmed, na.rm=T)

message2<-paste("The min is: ",min)
print (message2)

message3<-paste("The max is: ",max)
print (message3)

fivepercent = quantile(rowmed, probs = 0.03, na.rm=T)
message4<-paste("3% is: ",fivepercent)
print (message4)

tenpercent = quantile(rowmed, probs = 0.10, na.rm=T)
message5<-paste("10% is: ",tenpercent)
print (message5)

NinetyPercent = quantile(rowmed, probs = 0.90, na.rm=T)
message5<-paste("90% is: ",NinetyPercent)
print (message5)

NinetySevenPercent = quantile(rowmed, probs = 0.97, na.rm=T)
message5<-paste("97% is: ",NinetySevenPercent)
print (message5)
