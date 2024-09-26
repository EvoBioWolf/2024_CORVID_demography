
## parametric bootstrap with folded frequency spectrum

folded <- TRUE

source("jaatha_frsp_scrm_func.R")

i <- as.numeric(commandArgs()[6])

set.seed(123)

folded <- TRUE
origstats_folded <- unlist(read.table(paste0("fs_sim_folded/frsp_jss_folded_",i)))
jaatha_data_folded <- create_jaatha_data(origstats_folded, jm)                               
bsfrspf <- jaatha(jm, jaatha_data_folded, sim=1000, cores=8, repetitions=5)
save(bsfrspf, file=paste0("jfrsp_parbs_", i,"_123_folded.RData"))
