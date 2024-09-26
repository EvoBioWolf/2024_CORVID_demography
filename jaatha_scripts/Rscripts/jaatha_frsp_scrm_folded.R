folded <- TRUE
fullsize <- FALSE

source("jaatha_frsp_scrm_func.R")

## for testing: set.seed(123); estimatesfrsp <- jaatha(jm, jaatha_data, sim=40, cores=1, init_method = "middle")

set.seed(123)

estimatesfrsp <- jaatha(jm, jaatha_data_folded, sim=1000, cores=15, repetitions=5)

save(estimatesfrsp, file="jfrsp_zi_123_folded.RData")
