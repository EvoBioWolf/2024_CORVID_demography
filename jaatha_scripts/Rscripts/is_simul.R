library(scrm)
library(coala)

folded <- TRUE

source("jaatha_frsp_scrm_func.R")

indexrange <- 1:10 + 1000


load("jaatha_estimations/jfrsp_zi_123_folded.RData")
frsp_simres <- list()
for(i in indexrange) {
    frsp_simres[[i]] <- data_simulator(estimatesfrsp$estimate, fullsize=TRUE, return_also_jsfs=TRUE)
    cat(frsp_simres[[i]][[1]], "\n", file=paste0("is_sim_folded/frsp_jss_folded_", i))
    cat(frsp_simres[[i]][[2]], sep="\n", file=paste0("is_sim_folded/frsp_simspec_folded_", i))
}
