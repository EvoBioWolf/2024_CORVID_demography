library(usethis)
library(devtools) 
library(igraph)
library(admixtools)
library(tidyverse)
library(magrittr)
library(dplyr)
library(lobstr, lib="/dss/dsshome1/lxc0E/di67kah/R")
library(foreach, lib="/dss/dsshome1/lxc0E/di67kah/R")
library(iterators, lib = "/dss/dsshome1/lxc0E/di67kah/R")
library(parallel)
library(doParallel, lib="/dss/dsshome1/lxc0E/di67kah/R")
library(quadprog, lib="/dss/dsshome1/lxc0E/di67kah/R")
library(data.tree, lib="/dss/dsshome1/lxc0E/di67kah/R")
library(future, lib="/dss/dsshome1/lxc0E/di67kah/R")
library(furrr, lib="/dss/dsshome1/lxc0E/di67kah/R")

args <- commandArgs(trailingOnly = TRUE)
setwd("/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2/06_results/admixtools")

genotype_data = args[1]
dir = args[1]
# extract_f2(args[1], args[1])
mypops = c("am","cor1","cor2","cor3","cnx1","cnx2","cnx3P","cnx3S","cnx4","cnx5","cnx6","ori1","ori2","ori3","pec1")
f2_blocks = f2_from_precomp(dir, afprod = TRUE, pops=mypops)

df <- data.frame(run_no=c(1), bestscore=c(2))
 for (i in args[3]:args[4]) {
  opt_results = find_graphs(f2_blocks, outpop = "am", numadmix = as.numeric(args[2]),
                                stop_gen = 1000, stop_gen2=50) 
  winner = opt_results %>% slice_min(score, with_ties = FALSE)
  winner$score[[1]]
  df[i,1] <- i 
  df[i,2] <- winner$score[[1]]
  saveRDS(winner, file = paste(args[1],"_allpop_adm",args[2],"_run_",i, sep=""),
          ascii = FALSE, version = NULL,
          compress = TRUE, refhook = NULL)

  file_name = paste(args[1],"_allpop_adm",args[2],"_run_",i,".tiff", sep="")
  tiff(file_name, units="in", width=5, height=5, res=300)
  print(plot_graph(winner$edges[[1]]))
  dev.off()
  }
 print(df)
 write.table(df, file=paste(args[1],"_allpop_adm",args[2],"_bs",args[3],"to",args[4],"_bestscores.txt",sep=""), quote=FALSE, row.names = FALSE, sep="\t")
