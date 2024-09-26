
library(scrm)
library(coala)

source("../jaatha/sumstat_funcs.R")
source("param_conv.R")

scrm_frhd_simul_T <- function(par) {
    ## parameter conversion to ms or scrm in time scaling or 4NS generations,
    ## where NS is the current effective size of the spanish population; population sizes relative to NS
    ## see https://tskit.dev/msprime/docs/stable/switch_from_other_simulators.html
    
    L <- 10000
    NS <- par[11]
    fNS <- 4*NS
    
    par[1]/NS  ## Population size of ancestral population ASI of S and I.
    par[2]/fNS ## Time of the split of S and I.
    par[3]/NS  ## Population size of H right before H-W split
    par[4]/NS  ## Population size of H-I ancestral after this split but before H-I split
    par[5]/fNS ## Time T_H of split between H and I.
    par[6]/NS  ## Current size if H
    par[7] * fNS ## Growth rate of H since split (implicitly determines size of H at time T_H)
    par[8]/NS  ## Current size of I
    par[9]/fNS ## Time of the split of W from H
    par[10]/NS ## Current size of W
    par[12]*fNS  ## Migration rate from S to W (actually trace-back rate of lineages from W back to S)
    par[13]*fNS  ## Migration rate from W to S (actually trace-back rate of lineages from S back to W)
    par[14]*fNS  ## Migration rate from H to W (actually trace-back rate of lineages from W back to H)
    par[15]*fNS  ## Migration rate from W to H (actually trace-back rate of lineages from H back to W)
    par[16]*fNS  ## Migration rate from H to I (actually trace-back rate of lineages from I back to H)
    par[17]*fNS  ## Migration rate from I to H (actually trace-back rate of lineages from H back to I)
    par[18]/fNS  ## Time when migration between W and S started
    par[19] * fNS * (L-1)  ## recombination rate
    par[20]/fNS  ## time when growth of hooded crow population started
    par[21]/fNS  ## time when growth of W population started
    par[22]*fNS  ## growth rate of W
    3.18e-9 * fNS * L ## theta (per locus) for the Spanish population
    
    sampsize  <- 2*c(15, 24, 54, 5)  ## S, W, H, I
    
    ## growth rate needs extra treatment as scrm seems to have problems with growth rates close to 0
    gH <- par[7] * fNS 
    if(gH > -1e-5 & gH < 1e-5) gH <- 0
    gW <- par[22] * fNS 
    if(gW > -1e-5 & gW < 1e-5) gW <- 0
    
    s <- scrm(paste(sum(sampsize), 91,
                    "-T",
                    "-r", par[19] * fNS * (L-1), L,
                    "-I 4 30 48 108 10",
                    "-n 2", par[10]/NS,
                    "-n 3", par[6]/NS,
                    "-n 4", par[8]/NS,
                    "-g 3", gH,
                    "-g 2", gW,
                    "-m 1 2", par[13]*fNS,
                    "-m 2 1", par[12]*fNS,
                    "-m 2 3", par[14]*fNS,
                    "-m 3 2", par[15]*fNS,
                    "-m 3 4", par[17]*fNS,
                    "-m 4 3", par[16]*fNS,
                    "-em", par[18]/fNS, "1 2", 0,
                    "-em", par[18]/fNS, "2 1", 0,
                    "-ej", par[9]/fNS, "2 3",
                    "-em", par[9]/fNS, "2 3 0",
                    "-em", par[9]/fNS, "3 2 0",
                    "-en", par[9]/fNS, "3", par[3]/NS,
                    "-eg", par[20]/fNS, "3", 0,
                    "-eg", par[21]/fNS, "2", 0,
                    "-ej", par[5]/fNS, "4 3",
                    "-em", par[5]/fNS, "4 3 0",
                    "-em", par[5]/fNS, "3 4 0",
                    "-en", par[5]/fNS, "3", par[4]/NS,
                    "-ej", par[2]/fNS, "1 3",
                    "-em", par[2]/fNS, "1 3 0",
                    "-em", par[2]/fNS, "3 1 0",
                    "-en", par[2]/fNS, "3", par[1]/NS
                    ))
    s
}

scrm_frsp_simul_T <- function(par) {
    ## parameter conversion to ms or scrm in time scaling or 4NS generations,
    ## where NS is the current effective size of the spanish population; population sizes relative to NS
    ## see https://tskit.dev/msprime/docs/stable/switch_from_other_simulators.html
    
    L <- 10000
    NS <- par[11]
    fNS <- 4*NS
    
    par[1]/NS   ## Population size of ancestral population ASI of S and I.
    par[2]/fNS   ## Time of the split of S and I.
    par[3]/NS   ## Population size of S after this split but before S-W split
    par[4]/NS   ## Population size of I after this split but before H-I split
    par[5]/fNS   ## Time T_H of split between H and I.
    par[6]/NS   ## Current size if H
    par[7] * fNS  ## Growth rate of H since split (implicitly determines size of H at time T_H)
    par[8]/NS   ## Current size of I
    par[9]/fNS   ## Time T_W of the split of W from S
    par[10]/NS  ## Current size of W
    par[12]*fNS  ## Migration rate from S to W (actually trace-back rate of lineages from W back to S)
    par[13]*fNS  ## Migration rate from W to S (actually trace-back rate of lineages from S back to W)
    par[14]*fNS  ## Migration rate from H to W (actually trace-back rate of lineages from W back to H)
    par[15]*fNS  ## Migration rate from W to H (actually trace-back rate of lineages from H back to W)
    par[16]*fNS  ## Migration rate from H to I (actually trace-back rate of lineages from I back to H)
    par[17]*fNS  ## Migration rate from I to H (actually trace-back rate of lineages from H back to I)
    par[18]/fNS  ## Time when migration between W and H started
    par[19]*fNS*(L-1)  ## recombination rate
    par[20]/fNS  ## time when growth of hooded crow population started
    par[21]/fNS  ## time when growth of W population started
    par[22]*fNS  ## growth rate of W
    3.18e-9 * fNS * L ## theta (per locus) for the Spanish population
    
    sampsize  <- 2*c(15, 24, 54, 5)  ## S, W, H, I
    
    ## growth rate needs extra treatment as  scrm seems to have problems with growth rates close to 0
    gH <- par[7] * fNS 
    if(gH > -1e-5 & gH < 1e-5) gH <- 0
    gW <- par[22] * fNS 
    if(gW > -1e-5 & gW < 1e-5) gW <- 0
    
    s <- scrm(paste(sum(sampsize), 91,
                    "-T",
                    "-r", par[19] * fNS * (L-1) , L,
                    ## "-l 0",   ###  use this line for using SMC'
                    "-I 4 30 48 108 10",
                    "-n 2", par[10]/NS,
                    "-n 3", par[6]/NS,
                    "-n 4", par[8]/NS,
                    "-g 3", gH,
                    "-g 2", gW,
                    "-m 1 2", par[13]*fNS,
                    "-m 2 1", par[12]*fNS,
                    "-m 2 3", par[14]*fNS,
                    "-m 3 2", par[15]*fNS,
                    "-m 3 4", par[17]*fNS,
                    "-m 4 3", par[16]*fNS,
                    "-em", par[18]/fNS, "3 2", 0,
                    "-em", par[18]/fNS, "2 3", 0,
                    "-ej", par[9]/fNS, "2 1",
                    "-em", par[9]/fNS, "1 2", 0,
                    "-em", par[9]/fNS, "3 2", 0,
                    "-en", par[9]/fNS, "1", par[3]/NS,
                    "-ej", par[5]/fNS, "4 3",
                    "-em", par[5]/fNS, "3 4", 0,
                    "-en", par[5]/fNS, "3", par[4]/NS,
                    "-eg", par[20]/fNS, "3", 0,
                    "-eg", par[21]/fNS, "2", 0,
                    "-ej", par[2]/fNS, "1 3",
                    ##   "-eG", par[2]/fNS, 0,
                    ##   "-eM", par[2]/fNS, 0,
                    "-en", par[2]/fNS, "3", par[1]/NS
                    ))
    s
}

scrmouto01 <- function(d) {
    D <- array("-", dim=c(196, 10000))
    for(l in d[2:197]) {
        s <- strsplit(l, split=" +")[[1]]
        D[as.integer(s[1]),] <- strsplit(s[2], split="")[[1]]
    }
    bi <- which(apply(D, 2, function(x) length(unique(x)))==2)
    t(matrix(as.numeric(t(D[,bi]) != D[1, bi]), ncol=196))
}


sim_frsp_stats <- function(par, seed, only_jsfs=FALSE) {
    ### needs that global variable folded has value TRUE or FALSE
    set.seed(seed)
    print(conv_jds2mspr_sp(par))
    A <- array(NA, dim=c(91, 14))
    if(!only_jsfs) {
      s <- scrm_frsp_simul_T(conv_jds2mspr_sp(par))
      length(s[[1]])
      for(i in 1:91) {	
        length(s[[1]][[i]])
        cat(s[[1]][[i]], sep="\n", "\n", file="tmptree")
        d <- system(paste("seq-gen -l 10000 -s ", 3.18e-9 * 4 * exp(par[11]),
                          " -M HKY -t 2 -p ", length(s[[1]][[i]])," tmptree"),
                    intern=TRUE)
        D <- scrmouto01(d)
        sampsize  <- c(15, 24, 54, 5)
        y <- c(fgcv(D[1:(2*sampsize[1]),]),
               fgcv(D[(2*sampsize[1]+1):(2*sum(sampsize[1:2])),]),
               fgcv(D[(2*sum(sampsize[1:2])+1):(2*sum(sampsize[1:3])),]),
               fgcv(D[(2*sum(sampsize[1:3])+1):(2*sum(sampsize)),]),
               jfgcv(D[1:(2*sum(sampsize[1:2])),], sampsize[1]),
               jfgcv(D[(2*sampsize[1]+1):(2*sum(sampsize[1:3])),], sampsize[2]),
               jfgcv(D[(2*sum(sampsize[1:2])+1):(2*sum(sampsize)),], sampsize[3]))
        A[i, ] <- y
      }
    }
    B <- A[,(1:7)*2-1]/(A[,(1:7)*2]*(A[,(1:7)*2]-1)/2)
    simst <- categorize_fgc(B)
    d <- system(paste("../msprime/frsp_genome.py",
                      paste(conv_jds2mspr_sp(par), collapse=" "),
                      seed, folded
                      ), intern = TRUE)
    if(folded) {
        i <- which(d=="joint allele frequency spectrum for populations S and W:")
        simst <- c(simst, coarse_jsfs(t(sapply(strsplit(d[(i+1):(i+31)], " "), as.numeric)), folded=TRUE))
        i <- which(d=="joint allele frequency spectrum for populations S and H:")
        simst <- c(simst, coarse_jsfs(t(sapply(strsplit(d[(i+1):(i+31)], " "), as.numeric)), folded=TRUE))
        i <- which(d=="joint allele frequency spectrum for populations S and I:")
        simst <- c(simst, coarse_jsfs(t(sapply(strsplit(d[(i+1):(i+31)], " "), as.numeric)), folded=TRUE))
        i <- which(d=="joint allele frequency spectrum for populations W and H:")
        simst <- c(simst, coarse_jsfs(t(sapply(strsplit(d[(i+1):(i+31)], " "), as.numeric)), folded=TRUE))
        i <- which(d=="joint allele frequency spectrum for populations W and I:")
        simst <- c(simst, coarse_jsfs(t(sapply(strsplit(d[(i+1):(i+31)], " "), as.numeric)), folded=TRUE))
        i <- which(d=="joint allele frequency spectrum for populations H and I:")
        simst <- c(simst, coarse_jsfs(t(sapply(strsplit(d[(i+1):(i+31)], " "), as.numeric)), folded=TRUE))
    } else {
        
        j <- which(d=="joint allele frequency spectrum for populations S and W:")
        simst <- c(simst, coarse_jsfs(t(sapply(strsplit(d[(j+1):(j+31)], " "), as.numeric)), cap=confuse_ancestral_prob))
        j <- which(d=="joint allele frequency spectrum for populations S and H:")
        simst <- c(simst, coarse_jsfs(t(sapply(strsplit(d[(j+1):(j+31)], " "), as.numeric)), cap=confuse_ancestral_prob))
        j <- which(d=="joint allele frequency spectrum for populations S and I:")
        simst <- c(simst, coarse_jsfs(t(sapply(strsplit(d[(j+1):(j+31)], " "), as.numeric)), cap=confuse_ancestral_prob))
        j <- which(d=="joint allele frequency spectrum for populations W and H:")
        simst <- c(simst, coarse_jsfs(t(sapply(strsplit(d[(j+1):(j+31)], " "), as.numeric)), cap=confuse_ancestral_prob))
        j <- which(d=="joint allele frequency spectrum for populations W and I:")
        simst <- c(simst, coarse_jsfs(t(sapply(strsplit(d[(j+1):(j+31)], " "), as.numeric)), cap=confuse_ancestral_prob))
        j <- which(d=="joint allele frequency spectrum for populations H and I:")
        simst <- c(simst, coarse_jsfs(t(sapply(strsplit(d[(j+1):(j+31)], " "), as.numeric)), cap=confuse_ancestral_prob))
    }
    list(simst, d)
}

sim_frhd_stats <- function(par, seed) {
    ### needs that global variable folded has value TRUE or FALSE
    set.seed(seed)
    print(conv_jds2mspr_hd(par))
    s <- scrm_frhd_simul_T(conv_jds2mspr_hd(par))
    length(s[[1]])
    A <- array(NA, dim=c(91, 14))
    for(i in 1:91) {
        length(s[[1]][[i]])
        cat(s[[1]][[i]], sep="\n", "\n", file="tmptree")
        d <- system(paste("seq-gen -l 10000 -s ", 3.18e-9 * 4 * exp(par[11]),
                          " -M HKY -t 2 -p ", length(s[[1]][[i]])," tmptree"),
                    intern=TRUE)
        D <- scrmouto01(d)
        sampsize  <- c(15, 24, 54, 5)
        y <- c(fgcv(D[1:(2*sampsize[1]),]),
               fgcv(D[(2*sampsize[1]+1):(2*sum(sampsize[1:2])),]),
               fgcv(D[(2*sum(sampsize[1:2])+1):(2*sum(sampsize[1:3])),]),
               fgcv(D[(2*sum(sampsize[1:3])+1):(2*sum(sampsize)),]),
               jfgcv(D[1:(2*sum(sampsize[1:2])),], sampsize[1]),
               jfgcv(D[(2*sampsize[1]+1):(2*sum(sampsize[1:3])),], sampsize[2]),
               jfgcv(D[(2*sum(sampsize[1:2])+1):(2*sum(sampsize)),], sampsize[3]))
        A[i, ] <- y
    }
    B <- A[,(1:7)*2-1]/(A[,(1:7)*2]*(A[,(1:7)*2]-1)/2)
    simst <- categorize_fgc(B)
    d <- system(paste("../msprime/frhd_genome.py",
                      paste(conv_jds2mspr_hd(par), collapse=" "),
                      seed, folded
                      ), intern = TRUE)
    if(folded) {
        i <- which(d=="joint allele frequency spectrum for populations S and W:")
        simst <- c(simst, coarse_jsfs(t(sapply(strsplit(d[(i+1):(i+31)], " "), as.numeric)), folded=TRUE))
        i <- which(d=="joint allele frequency spectrum for populations S and H:")
        simst <- c(simst, coarse_jsfs(t(sapply(strsplit(d[(i+1):(i+31)], " "), as.numeric)), folded=TRUE))
        i <- which(d=="joint allele frequency spectrum for populations W and I:")
        simst <- c(simst, coarse_jsfs(t(sapply(strsplit(d[(i+1):(i+31)], " "), as.numeric)), folded=TRUE))
        i <- which(d=="joint allele frequency spectrum for populations W and H:")
        simst <- c(simst, coarse_jsfs(t(sapply(strsplit(d[(i+1):(i+31)], " "), as.numeric)), folded=TRUE))
        i <- which(d=="joint allele frequency spectrum for populations W and I:")
        simst <- c(simst, coarse_jsfs(t(sapply(strsplit(d[(i+1):(i+31)], " "), as.numeric)), folded=TRUE))
        i <- which(d=="joint allele frequency spectrum for populations H and I:")
        simst <- c(simst, coarse_jsfs(t(sapply(strsplit(d[(i+1):(i+31)], " "), as.numeric)), folded=TRUE))
    } else {
        
        j <- which(d=="joint allele frequency spectrum for populations S and W:")
        simst <- c(simst, coarse_jsfs(t(sapply(strsplit(d[(j+1):(j+31)], " "), as.numeric)), cap=confuse_ancestral_prob))
        j <- which(d=="joint allele frequency spectrum for populations S and H:")
        simst <- c(simst, coarse_jsfs(t(sapply(strsplit(d[(j+1):(j+31)], " "), as.numeric)), cap=confuse_ancestral_prob))
        j <- which(d=="joint allele frequency spectrum for populations S and I:")
        simst <- c(simst, coarse_jsfs(t(sapply(strsplit(d[(j+1):(j+31)], " "), as.numeric)), cap=confuse_ancestral_prob))
        j <- which(d=="joint allele frequency spectrum for populations W and H:")
        simst <- c(simst, coarse_jsfs(t(sapply(strsplit(d[(j+1):(j+31)], " "), as.numeric)), cap=confuse_ancestral_prob))
        j <- which(d=="joint allele frequency spectrum for populations W and I:")
        simst <- c(simst, coarse_jsfs(t(sapply(strsplit(d[(j+1):(j+31)], " "), as.numeric)), cap=confuse_ancestral_prob))
        j <- which(d=="joint allele frequency spectrum for populations H and I:")
        simst <- c(simst, coarse_jsfs(t(sapply(strsplit(d[(j+1):(j+31)], " "), as.numeric)), cap=confuse_ancestral_prob))
    }
    list(simst, d)
}


indexrange  <- 1:100

folded <- TRUE

## CHANGE FALSE TO TRUE FOR THE ANALYSIS TO BE DONE


if(FALSE) {
    load("jaatha_estimations/jfrsp_zi_123_folded.RData")
    frsp_simres <- list()
    for(i in indexrange) {
        frsp_simres[[i]] <- sim_frsp_stats(estimatesfrsp$estimate, i)
        cat(frsp_simres[[i]][[1]], "\n", file=paste0("fs_sim_folded/frsp_jss_folded_", i))
        cat(frsp_simres[[i]][[2]], sep="\n", file=paste0("fs_sim_folded/frsp_simspec_folded_", i))
    }   
}

if(FALSE) {
    load("jaatha_estimations/jfrsp_zi_123_folded_bs_corrected.RData")
    frsp_simres <- list()
    for(i in indexrange) {
        frsp_simres[[i]] <- sim_frsp_stats(bscorf_raw, i)
        cat(frsp_simres[[i]][[1]], "\n", file=paste0("fs_sim_folded_bscor/frsp_jss_folded_bscor_", i))
        cat(frsp_simres[[i]][[2]], sep="\n", file=paste0("fs_sim_folded_bscor/frsp_simspec_folded_bscor_", i))
    }   
}


if(FALSE) {
load("jaatha_estimations/jfrhd_zi_123_folded.RData")
frhd_simres <- list()
for(i in indexrange) {
    frhd_simres[[i]] <- sim_frhd_stats(estimatesfrhd$estimate, i)
    cat(frhd_simres[[i]][[1]], "\n", file=paste0("fs_sim_folded/frhd_jss_folded_", i))
    cat(frhd_simres[[i]][[2]], sep="\n", file=paste0("fs_sim_folded/frhd_simspec_folded_", i))
}
}




