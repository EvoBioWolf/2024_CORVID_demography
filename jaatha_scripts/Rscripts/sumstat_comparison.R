
folded <- TRUE
source("calc_stats.R")
origstats

m4xstats <- array(dim=c(100,100), NA)
for(i in c(1:20, 31:50, 61:100)) {
    m4xstats[i,] <- unlist(as.vector(read.table(paste0("fs_sim_folded/m4x_jss_folded_", i))))
}

apply(m4xstats, 2, mean, na.rm=TRUE)

spstats <- array(dim=c(100,100), NA)
for(i in c(1:100)) {
    spstats[i,] <- unlist(as.vector(read.table(paste0("fs_sim_folded/frsp_jss_folded_", i))))
}

apply(spstats, 2, mean, na.rm=TRUE)

cbind(apply(m4xstats, 2, mean, na.rm=TRUE), origstats, apply(spstats, 2, mean, na.rm=TRUE))
