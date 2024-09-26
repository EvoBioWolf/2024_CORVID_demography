### calculations of the summary statistics for the real data

source("../Rscripts/sumstat_funcs.R")

vcf  <- readLines("../linkage_info/intergenic_block_10000_variantsitesonly.recode.vcf") ### REPLACE BY DATA INPUT FILE FOR LINKAGE BLOCKS

## read file about which sample is which population
popoft <- read.table("../jSFS/cor1_cor2to3_cnx1to3_cnx6.pop", stringsAsFactors = FALSE) ### REPLACE BY DATA INPUT FILE FOR SAMPLE INFORMATION
popof <- character()
popof[sapply(strsplit(popoft$V1, '_'), '[', 1)]  <- popoft$V2

l  <- list()
snppos <- list()
sampname <- strsplit(vcf[[1]],'\t')[[1]][10:107]
nsnps  <- table(sapply(vcf[2:length(vcf)], function(x) strsplit(x, '\t')[[1]][1]))
## number of SNPS per scaffold
for(sc in names(nsnps)) {
    l[[sc]] <- matrix(NA, ncol=nsnps[sc], nrow=length(sampname)*2)
}

lc <- integer()
lc[names(nsnps)] <- 1

for(i in 2:length(vcf)) {
    spl <- strsplit(vcf[[i]],'\t')[[1]]
    scn  <- spl[1]
    snppos[[scn]][lc[scn]] <- as.numeric(spl[2])
    l[[scn]][((1:98)-1)*2+1,lc[scn]] <- as.numeric(substr(spl[10:length(spl)],1,1))
    l[[scn]][(1:98)*2,lc[scn]] <- as.numeric(substr(spl[10:length(spl)],3,3))
    lc[scn] <- lc[scn] + 1
}

ssl <- list()

library(coala)

for(i in 1:length(snppos)) {
    snppos[[i]]  <- (snppos[[i]] - min(snppos[[i]])+0.5)/2000
    sel  <- apply(l[[i]], 2, any, na.rm=TRUE) 
    ssl[[i]] <- create_segsites(l[[i]][,sel], snppos[[i]][sel])
}

pop <- popof[rep(sampname, each=2)]

A <- array(NA, dim=c(91, 14))
for(i in 1:91) {
    y <- rep(NA, 14)
    y <- c(fgcv(ssl[[i]]$snps[pop=="cor1",]),
           fgcv(ssl[[i]]$snps[pop=="cor2to3",]),
           fgcv(ssl[[i]]$snps[pop=="cnx1to3",]),
           fgcv(ssl[[i]]$snps[pop=="cnx6",]),
           jfgcv(rbind(ssl[[i]]$snps[pop=="cor1", ], ssl[[i]]$snps[pop=="cor2to3",]),
                 sum(pop=="cor1")/2),
           jfgcv(rbind(ssl[[i]]$snps[pop=="cor2to3", ], ssl[[i]]$snps[pop=="cnx1to3",]),
                 sum(pop=="cor2to3")/2),
           jfgcv(rbind(ssl[[i]]$snps[pop=="cnx1to3", ], ssl[[i]]$snps[pop=="cnx6",]),
                 sum(pop=="cnx1to3")/2)
           )
    A[i, ] <- y
}

B <- A[,(1:7)*2-1]/(A[,(1:7)*2]*(A[,(1:7)*2]-1)/2)

origstats <- categorize_fgc(B)

## note that in 1.2 GB long genome, 16 Mio SNPs have been found and were reduced to approx 5 Mio SNPs

if(folded) {  ### REPLACE PATHS FOR INPUT FILES FOR PAIRWISE JOINT SITE-FREQUENCY SPECTRA
    LCS <- readLines("../jSFS/cor1_cor2to3_cnx1to3_cnx6_50ind_folded_all_nochr18/4PopModel1_jointMAFpop1_0.obs")
    ## apply t to make rows pop0 and cols pop1 and coarse the jsfs
    origstats <- c(origstats, coarse_jsfs(t(apply(sapply(strsplit(LCS, "[\t ]")[3:33], "[", 2:32), 1, as.numeric)), folded=TRUE))
    
    LCS <- readLines("../jSFS/cor1_cor2to3_cnx1to3_cnx6_50ind_folded_all_nochr18/4PopModel1_jointMAFpop2_0.obs")
    origstats <- c(origstats, coarse_jsfs(t(apply(sapply(strsplit(LCS, "[\t ]")[3:33], "[", 2:32), 1, as.numeric)), folded=TRUE))
    
    LCS <- readLines("../jSFS/cor1_cor2to3_cnx1to3_cnx6_50ind_folded_all_nochr18/4PopModel1_jointMAFpop3_0.obs")
    origstats <- c(origstats, coarse_jsfs(t(apply(sapply(strsplit(LCS, "[\t ]")[3:13], "[", 2:32), 1, as.numeric)), folded=TRUE))
    
    LCS <- readLines("../jSFS/cor1_cor2to3_cnx1to3_cnx6_50ind_folded_all_nochr18/4PopModel1_jointMAFpop2_1.obs")
    origstats <- c(origstats, coarse_jsfs(t(apply(sapply(strsplit(LCS, "[\t ]")[3:33], "[", 2:32), 1, as.numeric)), folded=TRUE))
    
    LCS <- readLines("../jSFS/cor1_cor2to3_cnx1to3_cnx6_50ind_folded_all_nochr18/4PopModel1_jointMAFpop3_1.obs")
    origstats <- c(origstats, coarse_jsfs(t(apply(sapply(strsplit(LCS, "[\t ]")[3:13], "[", 2:32), 1, as.numeric)), folded=TRUE))
    
    LCS <- readLines("../jSFS/cor1_cor2to3_cnx1to3_cnx6_50ind_folded_all_nochr18/4PopModel1_jointMAFpop3_2.obs")
    origstats <- c(origstats, coarse_jsfs(t(apply(sapply(strsplit(LCS, "[\t ]")[3:13], "[", 2:32), 1, as.numeric)), folded=TRUE))
    origstats_folded <- origstats
} else {
    LCS <- readLines("../jSFS/cor1_cor2to3_cnx1to3_cnx6_50ind_unfolded_all_nochr18/4PopModel1_jointDAFpop1_0.obs")
    ## apply t to make rows pop0 and cols pop1 and coarse the jsfs
    origstats <- c(origstats, coarse_jsfs(t(apply(sapply(strsplit(LCS, "[\t ]")[3:33], "[", 2:32), 1, as.numeric))))
    
    LCS <- readLines("../jSFS/cor1_cor2to3_cnx1to3_cnx6_50ind_unfolded_all_nochr18/4PopModel1_jointDAFpop2_0.obs")
    origstats <- c(origstats, coarse_jsfs(t(apply(sapply(strsplit(LCS, "[\t ]")[3:33], "[", 2:32), 1, as.numeric))))
    
    LCS <- readLines("../jSFS/cor1_cor2to3_cnx1to3_cnx6_50ind_unfolded_all_nochr18/4PopModel1_jointDAFpop3_0.obs")
    origstats <- c(origstats, coarse_jsfs(t(apply(sapply(strsplit(LCS, "[\t ]")[3:13], "[", 2:32), 1, as.numeric))))
    
    LCS <- readLines("../jSFS/cor1_cor2to3_cnx1to3_cnx6_50ind_unfolded_all_nochr18/4PopModel1_jointDAFpop2_1.obs")
    origstats <- c(origstats, coarse_jsfs(t(apply(sapply(strsplit(LCS, "[\t ]")[3:33], "[", 2:32), 1, as.numeric))))
    
    LCS <- readLines("../jSFS/cor1_cor2to3_cnx1to3_cnx6_50ind_unfolded_all_nochr18/4PopModel1_jointDAFpop3_1.obs")
    origstats <- c(origstats, coarse_jsfs(t(apply(sapply(strsplit(LCS, "[\t ]")[3:13], "[", 2:32), 1, as.numeric))))
    
    LCS <- readLines("../jSFS/cor1_cor2to3_cnx1to3_cnx6_50ind_unfolded_all_nochr18/4PopModel1_jointDAFpop3_2.obs")
    origstats <- c(origstats, coarse_jsfs(t(apply(sapply(strsplit(LCS, "[\t ]")[3:13], "[", 2:32), 1, as.numeric))))
    origstats_unfolded <- origstats
}
