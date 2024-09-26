library(jaatha)
library(scrm)

folded <- TRUE

source("calc_stats.R")
source("param_conv.R")
origstats

scrm_frhd_simul <- function(par) {
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
    
    sampsize  <- 2*c(15, 24, 54, 5)  ## S, C, H, I

    ## growth rate needs extra treatment as there sees to be a bug in scrm with growth rates close to 0
    gH <- par[7] * fNS 
    if(gH > -1e-5 & gH < 1e-5) gH <- 0
    gW <- par[22] * fNS 
    if(gW > -1e-5 & gW < 1e-5) gW <- 0
    
    s <- scrm(paste(sum(sampsize), 91,
                    "-t", 3.18e-9 * fNS * L,
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
                    "-en", par[9]/fNS, "3", par[3]/NS,
                    "-eg", par[20]/fNS, "3", 0,
                    "-eg", par[21]/fNS, "2", 0,
                    "-ej", par[5]/fNS, "4 3",
                    "-en", par[5]/fNS, "3", par[4]/NS,
                    "-ej", par[2]/fNS, "1 3",
                    "-en", par[2]/fNS, "3", par[1]/NS
                    ))
    s
}

### try out: scrm_frhd_simul(c(166816.575522853, 153855.36170506, 3527.91059843214, 473080.71841167, 77684.4848921399, 10284.9360647968, 0.00388917469680905, 3964.74648204392, 73616.0157580889, 7866.94886883787, 3584.11048259564, 0.000592383793332374, 0.000248982487338623, 0.000116104256314031, 0.000219972362373255, 9.1094254550228e-05, 0.000682478109993667, 424.218735072772, 4.40701584569285e-10, 2875955316))
## 
### scrm_frhd_simul(c(10.17, 7.515, 9.554, 10.17, 8.204, 10.17, 0.0009977, 10.17, 8.204, 9.554, 9.554, 0.00045, 0.00045, 0.00045, 0.00045, 0.00045, 0.00045, 6.586, -19.989, 123))

data_simulator <- function(par, fullsize=FALSE) {
    ## global variable folded must be TRUE or FALSE
    ## entries of par are:
    ## 1: log of Population size of ancestral population SI of S and I.
    ## 2: log of Time between the S--I split and the H--I split.
    ## 3: log of size of HW at the time of the H--W split
    ## 4: log of Population size of HI, that is I after this split but before H-I split
    ## 5: log of Time between H--I and H--W split
    ## 6: log of Population size of H at present time
    ## 7: Growth rate of H (and HI) since H--I split
    ## 8: log of Current size of I
    ## 9: log of Time betweeen the split of W from H and the time when migration between W and S started
    ## 10: log of Current size of W
    ## 11: log of Current size of S
    ## 12: Migration rate from S to W (actually trace-back rate of lineages from W back to S)
    ## 13: Migration rate from W to S (actually trace-back rate of lineages from S back to W)
    ## 14: Migration rate from H to W (actually trace-back rate of lineages from W back to H)
    ## 15: Migration rate from W to H (actually trace-back rate of lineages from H back to W)
    ## 16: Migration rate from H to I (actually trace-back rate of lineages from I back to H)
    ## 17: Migration rate from I to H (actually trace-back rate of lineages from H back to I)
    ## 18: log of the time when the migration between W and S started
    ## 19: log of recombination rate
    ## 20: fraction of time of growth of hooded crow population compared to existence since split from W
    ## 21: fraction of time of growth of W population relative to existence since split from S
    ## 22: growth rate of W
    
    ## folded indicates whether the JSFSs should be folded
    
##    filename <- paste(c("frhd_simruns/scrm",signif(par, 4)), collapse="_")
    
##    system(paste("touch", filename))
    ## start_time <- Sys.time()
    s <- scrm_frhd_simul(conv_jds2mspr_hd(par))

    ## end_time  <- Sys.time()
##    system(paste("rm", filename))
    ## print(end_time - start_time)
    
    sampsize  <-  c(15, 24, 54, 5)
    nseq  <- 2*sum(sampsize)
    A <- array(NA, dim=c(91, 14))
    for(i in 1:91) {
        D <- s[[1]][[i]]
        if(ncol(D) < 2) {
            y <- rep(0,14)
        } else {
            y <- c(fgcv(D[1:(2*sampsize[1]),]),
                   fgcv(D[(2*sampsize[1]+1):(2*sum(sampsize[1:2])),]),
                   fgcv(D[(2*sum(sampsize[1:2])+1):(2*sum(sampsize[1:3])),]),
                   fgcv(D[(2*sum(sampsize[1:3])+1):(2*sum(sampsize)),]),
                   jfgcv(D[1:(2*sum(sampsize[1:2])),], sampsize[1]),
                   jfgcv(D[(2*sampsize[1]+1):(2*sum(sampsize[1:3])),], sampsize[2]),
                   jfgcv(D[(2*sum(sampsize[1:2])+1):(2*sum(sampsize)),], sampsize[3]))}
        A[i,] <- y
    }
    B <- A[,1:7*2-1]/(A[,1:7*2]*(A[,1:7*2]-1)/2)
    simst <- categorize_fgc(B)

##    filename <- paste(c("frhd_simruns/msprime",signif(par, 4)), collapse="_")
    
##    system(paste("touch", filename))
    ##start_time <- Sys.time()
    d <- system(paste("../msprime/fromhooded_no_rcmb.py",
                      paste(conv_jds2mspr_hd(par), collapse=" "),
                      sprintf("%12.0f", round(runif(1, 0, 4e9))), fullsize, folded
                      ), intern = TRUE)
    ##end_time <- Sys.time()
##    system(paste("rm", filename))
    ##cat(c(par, end_time - start_time, "\n"), file="frhd_runtimes", append=TRUE)  

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
        if(!fullsize) simst[24:length(simst)] <-  simst[24:length(simst)] * (625134484+11146221)/8e5
        ## corrected scaling as scaling below did not take into account that 14.823e6 SNPs are already without those in repeat regions.
        ## if(!fullsize) simst[24:length(simst)] <-  simst[24:length(simst)] * 1.2e9*11146221/14.823e6/8e5
        ## scaling: total genome length 1.2 GB, approx 11,146,221 of 14.823e6 SNPs have been chosen, 400*2000=8e5 positions have been simulated
        ## Check file WHATIS
        
    } else {
        
        i <- which(d=="joint allele frequency spectrum for populations S and W:")
        simst <- c(simst, coarse_jsfs(t(sapply(strsplit(d[(i+1):(i+31)], " "), as.numeric)), cap=confuse_ancestral_prob))
        i <- which(d=="joint allele frequency spectrum for populations S and H:")
        simst <- c(simst, coarse_jsfs(t(sapply(strsplit(d[(i+1):(i+31)], " "), as.numeric)), cap=confuse_ancestral_prob))
        i <- which(d=="joint allele frequency spectrum for populations S and I:")
        simst <- c(simst, coarse_jsfs(t(sapply(strsplit(d[(i+1):(i+31)], " "), as.numeric)), cap=confuse_ancestral_prob))
        i <- which(d=="joint allele frequency spectrum for populations W and H:")
        simst <- c(simst, coarse_jsfs(t(sapply(strsplit(d[(i+1):(i+31)], " "), as.numeric)), cap=confuse_ancestral_prob))
        i <- which(d=="joint allele frequency spectrum for populations W and I:")
        simst <- c(simst, coarse_jsfs(t(sapply(strsplit(d[(i+1):(i+31)], " "), as.numeric)), cap=confuse_ancestral_prob))
        i <- which(d=="joint allele frequency spectrum for populations H and I:")
        simst <- c(simst, coarse_jsfs(t(sapply(strsplit(d[(i+1):(i+31)], " "), as.numeric)), cap=confuse_ancestral_prob))
        if(!fullsize) simst[24:length(simst)] <-  simst[24:length(simst)] * (625134484+11146221)/8e5
        ## corrected scaling as scaling below did not take into account that 14.823e6 SNPs are already without those in repeat regions.
        ## if(!fullsize) simst[24:length(simst)] <-  simst[24:length(simst)] * 1.2e9*9265471/14.823e6/8e5
        ## scaling: total genome length 1.2 GB, approx 9,265,471 of 14.823e6 SNPs have been chosen, 400*2000=8e5 positions have been simulated
        ## Check file WHATIS 
    }
    simst
}

simulator <- function(par) {
    data_simulator(par)
}

id <- function(x,opts){x}

jm <- create_jaatha_model(simulator, par_ranges=matrix(c(
                                         log(1000), log(500000) , ## range lNeASI
                                         log(100), log(100000) ,  ## range lTsHISI"
                                         log(1000), log(500000) , ## range lNeAHC
                                         log(1000), log(500000) , ## range lNeAHI
                                         log(100), log(100000) ,  ## range lTsHCHI
                                         log(1000), log(500000) , ## range lNeH 
                                         0, 0.0001 ,            ## range gNeH 
                                         log(1000), log(500000) , ## range lNeI 
                                         log(100), log(100000) ,  ## range lTsHCCS
                                         log(1000), log(500000) , ## range lNeC 
                                         log(1000), log(500000) , ## range lNeS 
                                         0, 1e-3 , ## range "MSC"  
                                         0, 1e-3 , ## range "MCS"  
                                         0, 1e-3 , ## range "MHC"  
                                         0, 1e-3 , ## range "MCH"  
                                         0, 1e-3 , ## range "MHI"  
                                         0, 1e-3 , ## range "MIH"    
                                         log(100), log(10000),  ## range "lTmSC"  
                                         log(3.18e-9/3), log(3.18e-9*30), ## range "lrecr"
                                         0,1, ## range frTgH
                                         0,1, ## range frTgW
                                         0, 0.0001) , ## range "gNeW" 
                                         ncol=2, byrow=TRUE),
                          list(create_jaatha_stat("lNeASI", id),
                               create_jaatha_stat("lTsHISI" , id),
                               create_jaatha_stat("lNeAHC", id),
                               create_jaatha_stat("lNeAHI", id),
                               create_jaatha_stat("lTsHCHI" , id),
                               create_jaatha_stat("lNeH" , id),
                               create_jaatha_stat("gNeH" , id),
                               create_jaatha_stat("lNeI" , id),
                               create_jaatha_stat("lTsHCCS", id),
                               create_jaatha_stat("lNeC" , id),
                               create_jaatha_stat("lNeS" , id),
                               create_jaatha_stat("MSC"  , id),
                               create_jaatha_stat("MCS"  , id),
                               create_jaatha_stat("MHC"  , id),
                               create_jaatha_stat("MCH"  , id),
                               create_jaatha_stat("MHI"  , id),
                               create_jaatha_stat("MIH"  , id),
                               create_jaatha_stat("lTmSC"  , id),
                               create_jaatha_stat("lrecr", id),
                               create_jaatha_stat("frTgH", id),
                               create_jaatha_stat("frTgW", id),
                               create_jaatha_stat("gNeW", id)
                             )
                          )

if(folded) {
    jaatha_data_folded <- create_jaatha_data(origstats_folded, jm)
} else {
    jaatha_data_unfolded <- create_jaatha_data(origstats_unfolded, jm)
}
