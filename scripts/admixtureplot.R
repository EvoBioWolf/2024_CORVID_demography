library(sysfonts) #Rfont
library(showtextdb) #Rfont
library(showtext) #Rfont
par(family  = "Arial") #Rfont
showtext_auto() #Rfont

args <- commandArgs(trailingOnly = TRUE)
setwd(args[3])

safe_colorblind_palette <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", 
                             "#AA4499", "#44AA99", "#999933", "#882255", "#661100", 
                             "#6699CC", "#888888", "#000000", "#E69F00", "#56B4E9")

pdf(file=paste(args[1],"_Admixture_K=",args[2],".pdf",sep=""), width=8,height=4)
#(res=150, width=1400, height=600)
tbl <- read.table(file=paste(args[1],args[2],"Q",sep=".")) 
poplist <- read.table(args[4])
tbl$sample <- poplist$V4
tbl$pop <- poplist$V2
tbl$group <- poplist$V3
tbl.sorted <- tbl[order(tbl$group),]
barplot(t(as.matrix(tbl.sorted)), names.arg=tbl.sorted$sample, col=safe_colorblind_palette[1:10], xlab="", ylab="Ancestry", border=NA, las=2, cex.names=0.6)
dev.off()
