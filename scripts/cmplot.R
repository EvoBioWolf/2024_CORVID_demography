library(openxlsx)
library(CMplot) 
library(dplyr)
library(psych)
library(base)
library(sysfonts) #Rfont
library(showtextdb) #Rfont
library(showtext) #Rfont
par(family  = "Times") #Rfont
showtext_auto() #Rfont
options(digits = 9)
#conda R unable to find fonts in system. Using this as an alternative, otherwise redundant  

args <- commandArgs(trailingOnly = TRUE)

pat <- c("chr1","chr1A","chr2","chr3","chr4","chr4A","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12",
"chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chr23","chr24","chr25",
"chr26","chr27","chr28","chrZ","chrUN_1","chrUN_2","chrUN_3","chrUN_4","chrUN_5","chrUN_6","chrUN_7","chrUN_8","chrUN_9",
"chrUN_10","chrUN_11","chrUN_12","chrUN_13","chrUN_14","chrUN_15","chrUN_16","chrUN_17")

label <- c("1","1A","2","3","4","4A","5","6","7","8","9","10","11","12","13","14","15","17","18","19","20","21","22","23","24","26","27","28")

patno18 <- c("chr1","chr1A","chr2","chr3","chr4","chr4A","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12",
             "chr13","chr14","chr15","chr16","chr17","chr19","chr20","chr21","chr22","chr23","chr24","chr25",
             "chr26","chr27","chr28","chrZ","chrUN_1","chrUN_2","chrUN_3","chrUN_4","chrUN_5","chrUN_6","chrUN_7","chrUN_8","chrUN_9",
             "chrUN_10","chrUN_11","chrUN_12","chrUN_13","chrUN_14","chrUN_15","chrUN_16","chrUN_17")

scaffto18 <- c("scaffold_1026","scaffold_107","scaffold_78","scaffold_60","scaffold_750","scaffold_246","scaffold_927","scaffold_1223","scaffold_373","scaffold_70",
               "scaffold_261","scaffold_271","scaffold_305","scaffold_320",
               "scaffold_458","scaffold_971","scaffold_995","scaffold_1056")

setwd("/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2/06_results/fst")
for (i in args[1]) {
x <- read.delim(paste(args[2],"_",i,".win.chr", sep=""), header=FALSE, sep="")
x <- x[-1,]
dat <- data_frame() 
for (c in patno18) {
  chr <- paste("^",c,"$", sep="")
  mod <- x[grepl(chr, x$V2),]
  z <- seq(25000, length.out=nrow(mod), by=50000)
  modz <- cbind(z, mod)
  dat <- rbind(dat, modz)
}
mod <- x[grepl("^chr18$", x$V2),]
sorted.mod <- mod %>% arrange(factor(mod$V12, levels = scaffto18))
z <- seq(25000, length.out=nrow(mod), by=50000)
modz <- cbind(z, sorted.mod)
dat <- rbind(dat, modz)

sorted.x <- dat %>% arrange(factor(dat$V2, levels = pat))
sorted.x.2 <- sorted.x[!grepl("scaffold", sorted.x$V2),] #remove scaffolds not mapped to chr
sorted.x.chr <- sorted.x.2[!grepl("chrUN", sorted.x.2$V2),] #remove scaffolds not mapped to chr
fst <- as.character(sorted.x.chr$V10)
dxy <- as.character(sorted.x.chr$V9)
meanfst <- mean(as.numeric(fst), na.rm=TRUE)
print(meanfst)
thresdxy <- quantile(as.numeric(dxy), 0.999, na.rm=TRUE)
print(thresdxy)
thresfst <- quantile(as.numeric(fst), 0.999, na.rm=TRUE)
print(thresfst)

pi1 <- as.character(sorted.x.chr$V7)
pi2 <- as.character(sorted.x.chr$V8)
meanpi1 <- mean(as.numeric(pi1), na.rm=TRUE)
meanpi2 <- mean(as.numeric(pi2), na.rm=TRUE)

CMplot(sorted.x.chr[,c(2:3,1,10)],plot.type="m",LOG10=FALSE,col=c("grey30","grey60"),ylim=c(0,0.05),file="jpg",threshold=thresdxy,chr.labels.angle=0,
       signal.col=c("red"),signal.cex=c(0.5),memo=paste(i,"dxy_chr",sep="_"),dpi=150,file.output=TRUE,verbose=TRUE,cex=0.5,width=10,height=4,
       threshold.lty=2,threshold.lwd=1,threshold.col="red",ylab=paste("dxy: ",i,sep=""),cex.axis=1,chr.labels=label)

CMplot(sorted.x.chr[,c(2:3,1,12)],plot.type="m",LOG10=FALSE,col=c("grey30","grey60"),ylim=c(0,0.003),file="jpg",chr.labels.angle=0,
       signal.col=c("red"),signal.cex=c(0.5),memo=paste(i,"da_chr",sep="_"),dpi=150,file.output=TRUE,verbose=TRUE,cex=0.5,width=10,height=4,
       threshold.lty=2,threshold.lwd=1,ylab=paste("da: ",i, '\n', sep=""),cex.axis=1,chr.labels=label)

CMplot(sorted.x.chr[,c(2:3,1,11)],plot.type="m",LOG10=FALSE,col=c("grey30","grey60"),ylim=c(0,1),threshold=thresfst,file="jpg",chr.labels.angle=0,
       signal.col=c("red"),signal.cex=c(0.5),memo=paste(i,"fst_chr",sep="_"),dpi=150,file.output=TRUE,verbose=TRUE,cex=0.5,width=10,height=4,
       threshold.lty=2,threshold.lwd=1,threshold.col="red",ylab=paste("Fst: ",i, '\n', " (mean=", print(round(meanfst,4)), ")", sep=""),cex.axis=1,chr.labels=label)

CMplot(sorted.x.chr[,c(2:3,1,8:9)],plot.type="m",multracks=TRUE,LOG10=FALSE,col=c("grey60","#4197d8"),file="jpg",
       signal.col=NULL,signal.cex=c(1),memo=paste(i,"pi_chr",sep="_"),dpi=150,file.output=TRUE,verbose=TRUE,cex=0.5,width=10,height=4, highlight.text.cex=1.4,
       threshold.lty=2,threshold.lwd=1,ylab=paste("pi: ",i, '\n', " (mean=", print(round(meanpi1,4)), ", ",print(round(meanpi2,4)), ")", sep=""),chr.labels=label)


#specifically scaffold78: ~10M to 9M and scaffold60: 8M to 3.4M (aligning to ref CM02202.2), this rough offering is in reverse position and unknown scaffs are all dumped to the back (after scaff1026)

y <- read.delim(paste(args[2],"_",i,".win.chr18", sep=""), header=FALSE, sep="")
y <- y[-1,]
sorted.y <- y %>% arrange(factor(y$V12, levels = scaffto18))
z <- seq(25000, length.out=nrow(y), by=50000)
yz <- cbind(z, sorted.y)

write.table(yz, file=paste(args[2],"_",i, "_scafftochr18ref.txt", sep=""), sep = "\t", quote=FALSE,row.names = FALSE)

#defined elevated chr18 region as dummy pos 2525000 to 4275000
hd <- seq(2525000, 4275000 , by=50000)
yelevated <- data_frame() 
for (region in hd) 
{
  re <- paste("^",region,"$", sep="")
  elevated <- yz[grepl(re, yz$z),]
  yelevated <- rbind(yelevated, elevated)
}

write.table(yelevated, file=paste(args[2],"_", i, ".win.chr18.elevatedregion", sep=""),sep = "\t", quote=FALSE,row.names = FALSE)


CMplot(yz[,c(2:3,1,10)],plot.type="m",LOG10=FALSE,col=c("grey30","grey60"),ylim=c(0,0.02),file="jpg",
       signal.col=c("red"),signal.cex=c(1),memo=paste(i,"dxy_chr18",sep="_"),dpi=150,file.output=TRUE,verbose=TRUE,cex=0.5,width=10,height=4,
       ylab=paste("dxy: ",i,sep=""),cex.axis=1)

CMplot(yz[,c(2:3,1,12)],plot.type="m",LOG10=FALSE,col=c("grey30","grey60"),ylim=c(0,0.002),file="jpg",
       signal.col=c("red"),signal.cex=c(1),memo=paste(i,"da_chr18",sep="_"),dpi=150,file.output=TRUE,verbose=TRUE,cex=0.5,width=10,height=4,
       ylab=paste("da: ",i, '\n', sep=""),cex.axis=1)

CMplot(yz[,c(2:3,1,11)],plot.type="m",LOG10=FALSE,col=c("grey30","grey60"),ylim=c(0,1),file="jpg",
       signal.col=c("red"),signal.cex=c(1),memo=paste(i,"fst_chr18",sep="_"),dpi=150,file.output=TRUE,verbose=TRUE,cex=0.5,width=10,height=4,
       ylab=paste("Fst: ",i,sep=""),cex.axis=1)

CMplot(yz[,c(2:3,1,8:9)],plot.type="m",multracks=TRUE,LOG10=FALSE,ylim=c(0,0.02),col=c("grey60","#4197d8"),file="jpg",
       signal.col=NULL,signal.cex=c(1),dpi=150,memo=paste(i,"pi_chr18",sep="_"),file.output=TRUE,verbose=FALSE,cex=0.5,width=10,height=4,highlight.text.cex=1.4,
       threshold.lty=2,threshold.lwd=1,ylab=paste("pi: ",i, sep=""),chr.labels=label) 
}


