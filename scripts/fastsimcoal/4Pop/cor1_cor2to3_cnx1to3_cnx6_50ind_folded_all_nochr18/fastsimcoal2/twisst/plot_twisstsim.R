args <- commandArgs(trailingOnly = TRUE)

setwd("/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2/fastsimcoal2/4Pop/cor1_cor2to3_cnx1to3_cnx6_50ind_folded_all_nochr18/fastsimcoal2/bestruns/twisst")
setwd(paste(args[1],"tree_",args[2],sep=""))
source("/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2/06_results/twisst/plot_twisst.R")
weights_file <- paste(args[1],"_",args[2],".weights.csv.gz",sep="")
window_data_file <- paste(args[1],"_",args[2],".win.data.tsv",sep="")
twisst_data <- import.twisst(weights_files=weights_file, window_data_files=window_data_file)


############################## combined plots ##################################
#a summary plot shows all the topologies and a bar plot of their relative weightings
pdf(file=paste(args[1],"_",args[2],"_plotsummary.pdf",sep=""))
plot.twisst.summary(twisst_data, lwd=3, cex=1)
dev.off()
