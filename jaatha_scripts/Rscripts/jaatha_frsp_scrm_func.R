library(jaatha)
library(scrm)

source("../jaatha/calc_stats.R")
source("../jaatha/param_conv.R")
origstats

scrm_frsp_simul <- function(par) {
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
    par[7] * fNS  ## Growth rate of H (implicitly determines size of H at time when growth started)
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
    
    sampsize  <- 2*c(15, 24, 54, 5)  ## S, W, H, I   (W stands for Western Europe, including Germany and France)

    ## growth rate needs extra treatment as there sees to be a bug in scrm with growth rates close to 0
    gH <- par[7] * fNS
    if(gH > -1e-5 & gH < 1e-5) gH <- 0
    gW <- par[22] * fNS
    if(gW > -1e-5 & gW < 1e-5) gW <- 0
    
    
    s <- scrm(paste(sum(sampsize), 91,
                    "-t", 3.18e-9 * fNS * L,
                    "-r", par[19] * fNS * (L-1) , L,
                    ##   "-l 0",   ###  use this line for using SMC'
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
                    "-en", par[2]/fNS, "3", par[1]/NS
                    ))
    s
}

if(FALSE) {
    par <- c(10.015059,   8.059048,  10.015059,  10.015059,   8.059048,  10.015059,   0.000050,  10.015059,
             8.059048,  10.015059,  10.015059,   0.000500,   0.000500,   0.000500,   0.000500,   0.000500,
             0.000500,   6.907755, -18.415092,   0.500000,   0.500000,  0.000050)
    print(conv_jds2mspr_sp(par))
    for(i in 1:1000) {
        print(i)
        set.seed(i)
        start_time <- Sys.time()
        s <- scrm_frsp_simul(conv_jds2mspr_sp(par))
        end_time  <- Sys.time()
        print(end_time - start_time)
    }
}

data_simulator <- function(par, fullsize=FALSE, return_also_jsfs = FALSE) {
    ## note that global variable folded must be TRUE or FALSE
    
    ## entries of par are:
    ## 1: log of Population size of ancestral population ASI of S and I.
    ## 2: log of Time between the S--I split and the H--I split.
    ## 3: log of Population size of S after this split but before S-W split
    ## 4: log of Population size of I after this split but before H-I split
    ## 5: log of Time between H--I and S--W split
    ## 6: log of Current size if H
    ## 7: Growth rate of H
    ## 8: log of Current size of I (implicitly determins growth rate of I, assuming exponetial change)
    ## 9: log of Time betweeen the split of W from S and the time when migration between W and H started
    ## 10: log of Current size of W
    ## 11: log of Current size of S
    ## 12: Migration rate from S to W (actually trace-back rate of lineages from W back to S)
    ## 13: Migration rate from W to S (actually trace-back rate of lineages from S back to W)
    ## 14: Migration rate from H to W (actually trace-back rate of lineages from W back to H)
    ## 15: Migration rate from W to H (actually trace-back rate of lineages from H back to W)
    ## 16: Migration rate from H to I (actually trace-back rate of lineages from I back to H)
    ## 17: Migration rate from I to H (actually trace-back rate of lineages from H back to I)
    ## 18: log of the time when the migration between W and H started
    ## 19: log of recombination rate
    ## 20: fraction of time of growth of hooded crow population compared to existence since split from I
    ## 21: fraction of time of growth of W population relative to existence since split from S
    ## 22: growth rate of W
    
    ## folded indicates whether the JSFSs should be folded
    
    ##filename <- paste(c("frsp_simruns/scrm",signif(par, 4)), collapse="_")
    
    ##system(paste("touch", filename))
    ##start_time <- Sys.time()
    
    s <- scrm_frsp_simul(conv_jds2mspr_sp(par))
    
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
    
    ##filename <- paste(c("frsp_simruns/msprime",signif(par, 4)), collapse="_")
    
    ##system(paste("touch", filename))
    ## start_time <- Sys.time()
    
    d <- system(paste("../msprime/fromspain_no_rcmb.py",
                      paste(conv_jds2mspr_sp(par), collapse=" "),
                      sprintf("%12.0f", round(runif(1, 0, 4e9))), fullsize, folded
                      ), intern = TRUE)
    ## end_time <- Sys.time()
    ##system(paste("rm", filename))
    ## print(end_time - start_time)

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
    if(return_also_jsfs) return(list(simst, d))
    simst
}

simulator <- function(par) {
    data_simulator(par)
}

id <- function(x,opts){x}

jm <- create_jaatha_model(simulator,
                          par_ranges=matrix(c(
                              log(1000), log(500000) , ## range "lNeASI
                              log(100), log(100000) , ## range "lTsHISI"
                              log(1000), log(500000) , ## range "lNeASW
                              log(1000), log(500000) , ## range "lNeAHI
                              log(100), log(100000) , ## range "lTsSWHI"
                              log(1000), log(500000) , ## range "lNeH" 
                              0, 0.0001 , ## range "gNeH" 
                              log(1000), log(500000) , ## range "lNeI" 
                              log(100), log(100000) , ## range "lTsHWWS"
                              log(1000), log(500000) , ## range "lNeW" 
                              log(1000), log(500000) , ## range "lNeS" 
                              0, 1e-3 , ## range "MSW"  
                              0, 1e-3 , ## range "MWS"  
                              0, 1e-3 , ## range "MHW"  
                              0, 1e-3 , ## range "MWH"  
                              0, 1e-3 , ## range "MHI"  
                              0, 1e-3 , ## range "MIH"    
                              log(100), log(10000),  ## range "lTmHW"  
                              log(3.18e-9/3), log(3.18e-9*30), ## range "lrecr"
                              0,1, ## range frTgH
                              0,1, ## range frTgW
                              0, 0.0001) , ## range "gNeW" 
                              ncol=2, byrow=TRUE),
                          list(create_jaatha_stat("lNeASI", id),
                               create_jaatha_stat("lTsHISI" , id),
                               create_jaatha_stat("lNeASW", id),
                               create_jaatha_stat("lNeAHI", id),
                               create_jaatha_stat("lTsSWHI" , id),
                               create_jaatha_stat("lNeH" , id),
                               create_jaatha_stat("gNeH" , id),
                               create_jaatha_stat("lNeI" , id),
                               create_jaatha_stat("lTsHWWS", id),
                               create_jaatha_stat("lNeW" , id),
                               create_jaatha_stat("lNeS" , id),
                               create_jaatha_stat("MSW"  , id),
                               create_jaatha_stat("MWS"  , id),
                               create_jaatha_stat("MHW"  , id),
                               create_jaatha_stat("MWH"  , id),
                               create_jaatha_stat("MHI"  , id),
                               create_jaatha_stat("MIH"  , id),
                               create_jaatha_stat("lTmHW"  , id),
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
