folded <- TRUE

source("jaatha_frhd_scrm_func.R")

set.seed(123)

folded <- TRUE

estimatesfrhd <- jaatha(jm, jaatha_data_folded, sim=1000, cores=15, repetitions=5)

save(estimatesfrhd, file="jfrhd_zi_123_folded.RData")
