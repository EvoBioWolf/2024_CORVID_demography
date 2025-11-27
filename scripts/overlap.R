library(dplyr)
library(ggpubr)
library(devtools)
library(grid)
library(ggplot2)
#library(ggvenn)
library(VennDiagram)
library(base)

args <- commandArgs(trailingOnly = TRUE)

setwd("/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2/05.1_recal/overlap")

gatk <- read.table("134inds_DP3GQ0Miss10fullinfoQ100_pos.txt", sep="")
sam <- read.table("samtools_DP3GQ0Miss10Q30_pos.txt", sep="")
angsd <- read.table("134inds_rescaled_angsdrecalq30mindepth400genodepth3_min121_pos.txt", sep="")
poslist<- list(GATK=gatk$V1, Samtools=sam$V1, ANGSD=angsd$V1)
#venn.diagram(x=poslist, filename="OverlapVariants.pdf", output=TRUE)

overlap <- calculate.overlap(poslist)
summary(overlap)
overlap_all <- data.frame(LocusName=overlap$a5)

write.table(overlap_all, file="overlapping_variants_pos.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
