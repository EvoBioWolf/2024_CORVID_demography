library(devtools)
library(ggplot2)
library(SNPRelate)
library(openxlsx)
library(dplyr)
library(psych)
library(ggpubr)

args <- commandArgs(trailingOnly = TRUE)

setwd("/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2")
poplist <- read.delim(file=paste(args[4],".txt",sep=""))

setwd(args[2])

# PCA with SNPRelate #
# mean FST with SNPRelate #
colorBlindBlack8  <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
safe_colorblind_palette <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888", "#000000", "#E69F00", "#56B4E9")

vcf.fn <- (file=paste(args[1],".vcf.gz",sep=""))
snpgdsVCF2GDS(vcf.fn,"vcf.gds", method="copy.num.of.ref")
snpgdsSummary("vcf.gds")
vcf.gdsfile <- snpgdsOpen("vcf.gds")

# Fst estimation
v <- snpgdsFst(vcf.gdsfile, population=poplist$pop, method="W&C84", autosome.only=FALSE)
v$Fst
v$MeanFst
summary(v$FstSNP)

v2 <- snpgdsFst(vcf.gdsfile, population=poplist$pop, method="W&H02", autosome.only=FALSE)
v2$Fst
v2$MeanFst
v2$Beta
summary(v2$FstSNP)

vcf.pca <- snpgdsPCA(vcf.gdsfile, autosome.only=FALSE, missing.rate=NaN) ##use autosome.only=false to include all loci
names(vcf.pca)
#variance proportion (%)
pc.percent <- vcf.pca$varprop*100
head(round(pc.percent, 2))
print(pc.percent)
tab <- data.frame(sample.id = vcf.pca$sample.id,
EV1 = vcf.pca$eigenvect[,1], # the first eigenvector
EV2 = vcf.pca$eigenvect[,2], # the second eigenvector
EV3 = vcf.pca$eigenvect[,3], 
EV4 = vcf.pca$eigenvect[,4], 
EV5 = vcf.pca$eigenvect[,5], 
EV6 = vcf.pca$eigenvect[,6], 
EV7 = vcf.pca$eigenvect[,7], 
EV8 = vcf.pca$eigenvect[,8], 
EV9 = vcf.pca$eigenvect[,9], 
EV10 = vcf.pca$eigenvect[,10], 
stringsAsFactors = FALSE)
#extract vectors out for external plotting
tab$POP <- poplist$pop

setwd("/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2")
setwd(args[3])
write.csv(tab, file=paste(args[1],"snprelate.txt",sep="_"))

pdf(file=paste(args[1],"pc1vs2_pca.pdf",sep="_"))
ggplot(tab, group=POP) +
  geom_point(aes(x=EV1, y=EV2, color=POP, shape=POP), size=4) +
  scale_shape_manual(name="Population", values=c(15, 16, 17, 18, 3, 4, 9, 10, 11, 12, 13, 5, 15, 16, 17)) +
  scale_color_manual(name="Population", values=safe_colorblind_palette) +
  labs(x=paste("PC1 ",round(pc.percent[1],2),"%",sep=""), y=paste("PC2 ",round(pc.percent[2],2),"%",sep="")) +
  guides(color=guide_legend("Population"),fill=guide_legend("Population")) +
  ggtitle(args[1]) +
  theme_bw(base_size=12) +
  theme(axis.line = element_line(colour = "black"),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.border = element_blank(),
  panel.background = element_blank(),
  axis.title.x = element_text(size=15),
  axis.title.y = element_text(size=15),
  axis.text.x=element_text (size=12),
  axis.text.y=element_text (size=12),
  legend.text=element_text(size=12),
  legend.title=element_text(size=12)) 
dev.off()

# geom_text(aes(x=EV1, y=EV2, label=poplist$sample_id), check_overlap = TRUE) +

pdf(file=paste(args[1],"pc2vs3_pca.pdf",sep="_"))
ggplot(tab, group=POP) +
  geom_point(aes(x=EV2, y=EV3, color=POP, shape=POP), size=4) +
  scale_shape_manual(name="Population", values=c(15, 16, 17, 18, 3, 4, 9, 10, 11, 12, 13, 5, 15, 16, 17)) +
  scale_color_manual(name="Population", values=safe_colorblind_palette) +
  labs(x=paste("PC2 ",round(pc.percent[2],2),"%",sep=""), y=paste("PC3 ",round(pc.percent[3],2),"%",sep="")) +
  guides(color=guide_legend("Population"),fill=guide_legend("Population")) +
  ggtitle(args[1]) +
  theme_bw(base_size=12) +
  theme(axis.line = element_line(colour = "black"),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.border = element_blank(),
  panel.background = element_blank(),
  axis.title.x = element_text(size=15),
  axis.title.y = element_text(size=15),
  axis.text.x=element_text (size=12),
  axis.text.y=element_text (size=12),
  legend.text=element_text(size=12),
  legend.title=element_text(size=12)) 
dev.off()

pdf(file=paste(args[1],"pc2vs4_pca.pdf",sep="_"))
ggplot(tab, group=POP) +
  geom_point(aes(x=EV2, y=EV4, color=POP, shape=POP), size=4) +
  scale_shape_manual(name="Population", values=c(15, 16, 17, 18, 3, 4, 9, 10, 11, 12, 13, 5, 15, 16, 17)) +
  scale_color_manual(name="Population", values=safe_colorblind_palette) +
  labs(x=paste("PC2 ",round(pc.percent[2],2),"%",sep=""), y=paste("PC4 ",round(pc.percent[4],2),"%",sep="")) +
  guides(color=guide_legend("Population"),fill=guide_legend("Population")) +
  ggtitle(args[1]) +
  theme_bw(base_size=12) +
  theme(axis.line = element_line(colour = "black"),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.border = element_blank(),
  panel.background = element_blank(),
  axis.title.x = element_text(size=15),
  axis.title.y = element_text(size=15),
  axis.text.x=element_text (size=12),
  axis.text.y=element_text (size=12),
  legend.text=element_text(size=12),
  legend.title=element_text(size=12))
dev.off()

pdf(file=paste(args[1],"pc2vs5_pca.pdf",sep="_"))
ggplot(tab, group=POP) +
  geom_point(aes(x=EV2, y=EV5, color=POP, shape=POP), size=4) +
  scale_shape_manual(name="Population", values=c(15, 16, 17, 18, 3, 4, 9, 10, 11, 12, 13, 5, 15, 16, 17)) +
  scale_color_manual(name="Population", values=safe_colorblind_palette) +
  labs(x=paste("PC2 ",round(pc.percent[2],2),"%",sep=""), y=paste("PC5 ",round(pc.percent[5],2),"%",sep="")) +
  guides(color=guide_legend("Population"),fill=guide_legend("Population")) +
  ggtitle(args[1]) +
  theme_bw(base_size=12) +
  theme(axis.line = element_line(colour = "black"),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.border = element_blank(),
  panel.background = element_blank(),
  axis.title.x = element_text(size=15),
  axis.title.y = element_text(size=15),
  axis.text.x=element_text (size=12),
  axis.text.y=element_text (size=12),
  legend.text=element_text(size=12),
  legend.title=element_text(size=12))
dev.off()

##Nomissingness
setwd("/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2")
setwd(args[2])
vcf.pca <- snpgdsPCA(vcf.gdsfile, autosome.only=FALSE, missing.rate=0) ##use autosome.only=false to include all loci
names(vcf.pca)
#variance proportion (%)
pc.percent <- vcf.pca$varprop*100
head(round(pc.percent, 2))
print(pc.percent)
tab <- data.frame(sample.id = vcf.pca$sample.id,
EV1 = vcf.pca$eigenvect[,1], # the first eigenvector
EV2 = vcf.pca$eigenvect[,2], # the second eigenvector
EV3 = vcf.pca$eigenvect[,3],
EV4 = vcf.pca$eigenvect[,4],
EV5 = vcf.pca$eigenvect[,5],
EV6 = vcf.pca$eigenvect[,6],
EV7 = vcf.pca$eigenvect[,7],
EV8 = vcf.pca$eigenvect[,8],
EV9 = vcf.pca$eigenvect[,9],
EV10 = vcf.pca$eigenvect[,10], 
stringsAsFactors = FALSE)
#extract vectors out for external plotting
tab$POP <- poplist$pop

setwd("/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2")
setwd(args[3])
write.csv(tab, file=paste(args[1],"nomissing_snprelate.txt",sep="_"))

pdf(file=paste(args[1],"nomissing_pc1vs2_pca.pdf",sep="_"))
ggplot(tab, group=POP) +
  geom_point(aes(x=EV1, y=EV2, color=POP, shape=POP), size=4) +
  scale_shape_manual(name="Population", values=c(15, 16, 17, 18, 3, 4, 9, 10, 11, 12, 13, 5, 15, 16, 17)) +
  scale_color_manual(name="Population", values=safe_colorblind_palette) +
  labs(x=paste("PC1 ",round(pc.percent[1],2),"%",sep=""), y=paste("PC2 ",round(pc.percent[2],2),"%",sep="")) +
  guides(color=guide_legend("Population"),fill=guide_legend("Population")) +
  ggtitle(args[1]) +
  theme_bw(base_size=12) +
  theme(axis.line = element_line(colour = "black"),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.border = element_blank(),
  panel.background = element_blank(),
  axis.title.x = element_text(size=15),
  axis.title.y = element_text(size=15),
  axis.text.x=element_text (size=12),
  axis.text.y=element_text (size=12),
  legend.text=element_text(size=12),
  legend.title=element_text(size=12)) 
dev.off()

pdf(file=paste(args[1],"nomissing_pc2vs3_pca.pdf",sep="_"))
ggplot(tab, group=POP) +
  geom_point(aes(x=EV2, y=EV3, color=POP, shape=POP), size=4) +
  scale_shape_manual(name="Population", values=c(15, 16, 17, 18, 3, 4, 9, 10, 11, 12, 13, 5, 15, 16, 17)) +
  scale_color_manual(name="Population", values=safe_colorblind_palette) +
  labs(x=paste("PC2 ",round(pc.percent[2],2),"%",sep=""), y=paste("PC3 ",round(pc.percent[3],2),"%",sep="")) +
  guides(color=guide_legend("Population"),fill=guide_legend("Population")) +
  ggtitle(args[1]) +
  theme_bw(base_size=12) +
  theme(axis.line = element_line(colour = "black"),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.border = element_blank(),
  panel.background = element_blank(),
  axis.title.x = element_text(size=15),
  axis.title.y = element_text(size=15),
  axis.text.x=element_text (size=12),
  axis.text.y=element_text (size=12),
  legend.text=element_text(size=12),
  legend.title=element_text(size=12))
dev.off()

##Relatedness analysis
#IBD using PLINK's method of moments
# setwd("/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2")
# setwd(args[2])
# ibd <- snpgdsIBDMoM(vcf.gdsfile, sample.id=poplist$sample_id,
#     maf=0.001, missing.rate=0.05, autosome.only=FALSE)
# ibd.coeff <- snpgdsIBDSelection(ibd)
# requires snp.id

#within-and-between family relationship inferred by the King-robust method
# ibd.robust <- snpgdsIBDKING(vcf.gdsfile, sample.id=poplist$sample_id,
#     family.id=poplist$pop)
# fam <- snpgdsIBDSelection(ibd.robust)

# setwd("/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2")
# setwd(args[3])
# write.csv(ibd.coeff, file=paste(args[1],"idb_mom.txt",sep="_"))
# write.csv(fam, file=paste(args[1],"idb_king.txt",sep="_"))

