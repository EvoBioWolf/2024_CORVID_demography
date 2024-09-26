source("param_conv.R")

L <- system("ls jaatha_estimations/jfrhd*", intern=TRUE)

for(l in L) {
    load(l)
    cat("\nFile Name:", l, "\n")
    cat("CompositeLogLikelihood:", estimatesfrhd$loglikelihood, "\n")
    cat("Jaatha scaled parameters:", p <- estimatesfrhd$estimate, "\n")
    cat("msprime scaled parameters:", p <- conv_jds2mspr_hd(p), "\n")
    cat("fastsimcoal scaled parameters:", p <- conv_mspr2fsc_hd(p), "\n")
    cat("fastsimcoal scaled parameters again, with names:\n")
    print(p)
}

data.frame(signif(p, 3))

L <- system("ls jaatha_estimations/jfrsp*", intern=TRUE)

for(l in L) {
    load(l)
    cat("\nFile Name:", l, "\n")
    cat("CompositeLogLikelihood:", estimatesfrsp$loglikelihood, "\n")
    cat("Jaatha scaled parameters:", p <- estimatesfrsp$estimate, "\n")
    cat("msprime scaled parameters:", p <- conv_jds2mspr_sp(p), "\n")
    cat("fastsimcoal scaled parameters:", p <- conv_mspr2fsc_sp(p), "\n")
    cat("fastsimcoal scaled parameters again, with names:\n")
    print(signif(p, 3))
}

data.frame(signif(p, 3))

L <- system("ls jaatha_estimations/jfrhd*", intern=TRUE)

for(l in L) {
    load(l)
    cat(l, estimatesfrhd$loglikelihood, "\n")
}

L <- system("ls jaatha_estimations/jfrsp*", intern=TRUE)

for(l in L) {
    load(l)
    cat(l, estimatesfrsp$loglikelihood, "\n")
}

mod4x <- c(4715, 19227, 337358, 37171, 1124533, 782806, 74018, 1522, 3.67104e-04, 5.44207e-04, 0.0013670, 1.64834e-04, 5.63613e-06, 3.88323e-05, -6.41772e-04, 72248, 24365, 1629)
names(mod4x)  <- c("NPOP1", "NPOP2", "NPOP3", "NPOP4", "NANC1", "NANC2", "NANC3", "TMRMG", "MIG01R", "MIG10R", "MIG12R", "MIG21R", "MIG23R", "MIG32R", "GR2", "TDIV3", "TDIV2", "TDIV1")
mod4x_j <- conv_mspr2jds_hd(conv_fsc2mspr_hd(c(mod4x, exp(-2.046097e+01))))


folded  <- TRUE
source("jaatha_frsp_scrm_func.R") ## SPAIN
jstats <- create_jaatha_data(origstats, jm)
load("jaatha_estimations/jfrsp_zi_123_folded.RData")
estimate_llh(jm, jstats,  estimatesfrsp$estimate)$value  ## 
load("jaatha_estimations/jfrsp_zi_123_folded.RData")
estimate_llh(jm, jstats,  estimatesfrsp$estimate)$value  ## jaatha parameter range problem
load("jaatha_estimations/jfrsp_zi_123.RData")
estimate_llh(jm, jstats,  estimatesfrsp$estimate)$value ## 
load("jaatha_estimations/jfrsp_zi_123.RData")
estimate_llh(jm, jstats,  estimatesfrsp$estimate)$value ## 
source("jaatha_frhd_scrm_func.R") ## HOODED
jstats <- create_jaatha_data(origstats, jm)
load("jaatha_estimations/jfrhd_zi_123_folded.RData") 
estimate_llh(jm, jstats,  estimatesfrhd$estimate)$value ##
load("jaatha_estimations/jfrhd_zi_123_folded.RData") 
estimate_llh(jm, jstats,  estimatesfrhd$estimate)$value ## 
load("jaatha_estimations/jfrhd_zi_123.RData")
estimate_llh(jm, jstats,  estimatesfrhd$estimate)$value ##
load("jaatha_estimations/jfrhd_zi_123.RData")
estimate_llh(jm, jstats,  estimatesfrhd$estimate)$value ## jaatha parameter range problem

## have to increase the parameter range for fastsimcoal estimation
jmx <- create_jaatha_model(simulator, par_ranges=matrix(c(
                                         log(1000), log(600000) , ## range lNeASI
                                         log(100), log(100000) ,  ## range lTsHISI"
                                         log(1000), log(600000) , ## range lNeAHC
                                         log(1000), log(600000) , ## range lNeAHI
                                         log(100), log(100000) ,  ## range lTsHCHI
                                         log(1000), log(600000) , ## range lNeH 
                                         0, 0.001 ,            ## range gNeH 
                                         log(1000), log(600000) , ## range lNeI 
                                         log(100), log(100000) ,  ## range lTsHCCS
                                         log(1000), log(600000) , ## range lNeC 
                                         log(1000), log(600000) , ## range lNeS 
                                         0, 1e-2 , ## range "MSC"  
                                         0, 1e-2 , ## range "MCS"  
                                         0, 1e-2 , ## range "MHC"  
                                         0, 1e-2 , ## range "MCH"  
                                         0, 1e-2 , ## range "MHI"  
                                         0, 1e-2 , ## range "MIH"    
                                         log(100), log(10000),  ## range "lTmSC"  
                                         log(3.18e-9/3), log(3.18e-9*30)), ## range "lrecr"
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
                               create_jaatha_stat("lrecr", id)
                             )
                          )

estimate_llh(jmx, jstats,  mod4x_j)$value ##  

folded  <- FALSE
confuse_ancestral_prob <- 0.001594
source("jaatha_frsp_scrm_func.R")
jstats <- create_jaatha_data(origstats, jm)
load("jaatha_estimations/jfrsp_zi_123_folded.RData") 
estimate_llh(jm, jstats,  estimatesfrsp$estimate)$value    ## 
load("jaatha_estimations/jfrsp_zi_123_folded.RData")
estimate_llh(jm, jstats,  estimatesfrsp$estimate)$value ## jaatha parameter range problem
load("jaatha_estimations/jfrsp_zi_123.RData")
estimate_llh(jm, jstats,  estimatesfrsp$estimate)$value ##  
load("jaatha_estimations/jfrsp_zi_123.RData")
estimate_llh(jm, jstats,  estimatesfrsp$estimate)$value ## 
source("jaatha_frhd_scrm_func.R")
jstats <- create_jaatha_data(origstats, jm)
load("jaatha_estimations/jfrhd_zi_123_folded.RData")
estimate_llh(jm, jstats,  estimatesfrhd$estimate)$value  ## 
load("jaatha_estimations/jfrhd_zi_123_folded.RData")
estimate_llh(jm, jstats,  estimatesfrhd$estimate)$value  ## 
load("jaatha_estimations/jfrhd_zi_123.RData")
estimate_llh(jm, jstats,  estimatesfrhd$estimate)$value  ## 
load("jaatha_estimations/jfrhd_zi_123.RData")
estimate_llh(jm, jstats,  estimatesfrhd$estimate)$value ## jaatha parameter range problem

## have to increase the parameter range for fastsimcoal estimation:
jmx <- create_jaatha_model(simulator, par_ranges=matrix(c(
                                         log(1000), log(600000) , ## range lNeASI
                                         log(100), log(100000) ,  ## range lTsHISI"
                                         log(1000), log(600000) , ## range lNeAHC
                                         log(1000), log(600000) , ## range lNeAHI
                                         log(100), log(100000) ,  ## range lTsHCHI
                                         log(1000), log(600000) , ## range lNeH 
                                         0, 0.001 ,            ## range gNeH 
                                         log(1000), log(600000) , ## range lNeI 
                                         log(100), log(100000) ,  ## range lTsHCCS
                                         log(1000), log(600000) , ## range lNeC 
                                         log(1000), log(600000) , ## range lNeS 
                                         0, 1e-2 , ## range "MSC"  
                                         0, 1e-2 , ## range "MCS"  
                                         0, 1e-2 , ## range "MHC"  
                                         0, 1e-2 , ## range "MCH"  
                                         0, 1e-2 , ## range "MHI"  
                                         0, 1e-2 , ## range "MIH"    
                                         log(100), log(10000),  ## range "lTmSC"  
                                         log(3.18e-9/3), log(3.18e-9*30)), ## range "lrecr"
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
                               create_jaatha_stat("lrecr", id)
                             )
                          )

estimate_llh(jmx, jstats,  mod4x_j)$value ## -167354122


### Estimation from the model in which the Central Europeans vom from hooded crows

## with folded frequency spectrum
load("jaatha_estimations/jfrhd_zi_123_folded.RData")

par <- estimatesfrhd$estimate
par <- c(exp(par[1]), exp(par[2])+exp(par[5])+exp(par[9])+exp(par[18]), exp(par[3]), exp(par[4]),
                      exp(par[5])+exp(par[9])+exp(par[18]), exp(par[6]),
                      par[7], exp(par[8]), exp(par[9])+exp(par[18]), exp(par[10]), exp(par[11]),
         par[12], par[13], par[14], par[15], par[16], par[17], exp(par[18]), exp(par[19]))

L <- 10000
NS <- par[11]
fNS <- 4*NS

estimates_hd_mod <- data.frame(row.names=c("Ne S", "Ne C", "Ne H", "Growth rate Ne H", "Ne I", "Ne Anc H-C", "Ne Anc H-I",
                                           "Ne Anc S-C-H-I", "t C-H", "t H-I", "t S-(C-H-I)", "t start mig C-S",
                                           "mig C <- S", "mig S <- C", "mig C <- H", "mig H <- C", "mig H <- I", "mig I <- H",
                                           "rho/theta"), best.frhd.jaatharun=rep(0.0, 19))



## all rates an times in units of 4 * NS generations, where NS is the effective size of the Spanish population
## and stored in dataframe the per generation rates and absolute popsizes

(estimates_hd_mod["Ne S", "best.frhd.jaatharun"] <- NS)
par[1]/NS  ## Population size of ancestral population ASI of S and I, relative to present spanish pop
(estimates_hd_mod["Ne Anc S-C-H-I", "best.frhd.jaatharun"] <- par[1])
par[2]/fNS ## Time of the split of S and I.
(estimates_hd_mod["t S-(C-H-I)", "best.frhd.jaatharun"] <- par[2])
par[3]/NS  ## Population size of H right before H-C split, relative to present spanish pop
(estimates_hd_mod["Ne Anc H-C", "best.frhd.jaatharun"] <- par[3])
par[4]/NS  ## Population size of H-I ancestral after this split but before H-I split, relative to present spanish pop
(estimates_hd_mod["Ne Anc H-I", "best.frhd.jaatharun"] <- par[4])
par[5]/fNS ## Time T_H of split between H and I.
(estimates_hd_mod["t H-I", "best.frhd.jaatharun"] <- par[5])
par[6]/NS  ## Current size if H, relative to present spanish pop
(estimates_hd_mod["Ne H", "best.frhd.jaatharun"] <- par[6])
par[7]*fNS ## Growth rate of H since split (implicitly determines size of H at time T_H)
(estimates_hd_mod["Growth rate Ne H", "best.frhd.jaatharun"] <- par[7])
par[8]/NS  ## Current size of I
(estimates_hd_mod["Ne I", "best.frhd.jaatharun"] <- par[8])
par[9]/fNS ## Time of the split of C from H
(estimates_hd_mod["t C-H", "best.frhd.jaatharun"] <- par[9])
par[10]/NS ## Current size of C
(estimates_hd_mod["Ne C", "best.frhd.jaatharun"] <- par[10])
par[12]*fNS  ## Migration rate from S to C (actually trace-back rate of lineages from C back to S)
(estimates_hd_mod["mig C <- S", "best.frhd.jaatharun"] <- par[12])
par[13]*fNS  ## Migration rate from C to S (actually trace-back rate of lineages from S back to C)
(estimates_hd_mod["mig S <- C", "best.frhd.jaatharun"] <- par[13])
par[14]*fNS  ## Migration rate from H to C (actually trace-back rate of lineages from C back to H)
(estimates_hd_mod["mig C <- H", "best.frhd.jaatharun"] <- par[14])
par[15]*fNS  ## Migration rate from C to H (actually trace-back rate of lineages from H back to C)
(estimates_hd_mod["mig H <- C", "best.frhd.jaatharun"] <- par[15])
par[16]*fNS  ## Migration rate from H to I (actually trace-back rate of lineages from I back to H)
(estimates_hd_mod["mig I <- H", "best.frhd.jaatharun"] <-par[16] )
par[17]*fNS  ## Migration rate from I to H (actually trace-back rate of lineages from H back to I)
(estimates_hd_mod["mig H <- I", "best.frhd.jaatharun"] <- par[17])
par[18]/fNS  ## Time when migration between C and S started
(estimates_hd_mod["t start mig C-S", "best.frhd.jaatharun"] <- par[18])
par[19] * fNS * (L-1)  ## recombination rate per locus of 10 kb
(estimates_hd_mod["rho/theta", "best.frhd.jaatharun"] <- par[19]/3.18e-9)
3.18e-9 * fNS * L ## theta (per locus) for the Spanish population


estimatesfrhd$loglikelihood

### estimations from the model in which the Central Europeans come from Spain:

## with folded frequency spectrum

estimates_sp_mod <- data.frame(row.names=c("Ne S", "Ne C", "Ne H", "Growth rate Ne H", "Ne I", "Ne Anc S-C", "Ne Anc H-I",
                                           "Ne Anc S-C-H-I", "t S-C", "t H-I", "t (S-C)-(H-I)", "t start mig C-H",
                                           "mig C <- S", "mig S <- C", "mig C <- H", "mig H <- C", "mig H <- I", "mig I <- H",
                                           "rho/theta"), best.frsp.jaatharun=rep(0.0, 19))

load("jaatha_estimations/jfrsp_zi_123_folded.RData")

par <- estimatesfrsp$estimate
par <- c(exp(par[1]), exp(par[2])+exp(par[5])+exp(par[9])+exp(par[18]), exp(par[3]),
                               exp(par[4]), exp(par[5])+exp(par[9])+exp(par[18]), exp(par[6]),
                               par[7], exp(par[8]), exp(par[9])+exp(par[18]), exp(par[10]), exp(par[11]),
         par[12], par[13], par[14], par[15], par[16], par[17], exp(par[18]), exp(par[19]))
 
    L <- 10000
    NS <- par[11]
fNS <- 4*NS


## all population sizes are relative to (that is multiples of) the Spanish population.
## Rates and times are in units of 4N generations, where N is the effective size of the Spanish population,
## and per locus of 10 kb

estimates_sp_mod["Ne S", "best.frsp.jaatharun"] <- NS   
par[1]/NS   ## Population size of ancestral population ASI of S and I.
(estimates_sp_mod["Ne Anc S-C-H-I", "best.frsp.jaatharun"] <- par[1])
par[2]/fNS   ## Time of the split of S and I.
(estimates_sp_mod["t (S-C)-(H-I)", "best.frsp.jaatharun"] <- par[2])
par[3]/NS   ## Population size of S after this split but before S-C split
(estimates_sp_mod["Ne Anc S-C", "best.frsp.jaatharun"] <- par[3])
par[4]/NS   ## Population size of I after this split but before H-I split
(estimates_sp_mod["Ne Anc H-I", "best.frsp.jaatharun"] <- par[4])
par[5]/fNS   ## Time T_H of split between H and I.
(estimates_sp_mod["t H-I", "best.frsp.jaatharun"] <- par[5])
par[6]/NS   ## Current size if H
(estimates_sp_mod["Ne H", "best.frsp.jaatharun"] <- par[6])
par[7]*fNS  ## Growth rate of H since split (implicitly determines size of H at time T_H)   TODO  : Check growth rate scaling !!!!!!!!!
(estimates_sp_mod["Growth rate Ne H", "best.frsp.jaatharun"] <- par[7])
par[8]/NS   ## Current size of I
(estimates_sp_mod["Ne I", "best.frsp.jaatharun"] <- par[8])
par[9]/fNS   ## Time T_C of the split of C from S
(estimates_sp_mod["t S-C", "best.frsp.jaatharun"] <- par[9])
par[10]/NS  ## Current size of C
(estimates_sp_mod["Ne C", "best.frsp.jaatharun"] <- par[10])
par[12]*fNS  ## Migration rate from S to C (actually trace-back rate of lineages from C back to S)
(estimates_sp_mod["mig C <- S", "best.frsp.jaatharun"] <- par[12])
par[13]*fNS  ## Migration rate from C to S (actually trace-back rate of lineages from S back to C)
(estimates_sp_mod["mig S <- C", "best.frsp.jaatharun"] <- par[13])
par[14]*fNS  ## Migration rate from H to C (actually trace-back rate of lineages from C back to H)
(estimates_sp_mod["mig C <- H", "best.frsp.jaatharun"] <- par[14])
par[15]*fNS  ## Migration rate from C to H (actually trace-back rate of lineages from H back to C)
(estimates_sp_mod["mig H <- C", "best.frsp.jaatharun"] <- par[15])
par[16]*fNS  ## Migration rate from H to I (actually trace-back rate of lineages from I back to H)
(estimates_sp_mod["mig I <- H", "best.frsp.jaatharun"] <- par[16])
par[17]*fNS  ## Migration rate from I to H (actually trace-back rate of lineages from H back to I)
(estimates_sp_mod["mig H <- I", "best.frsp.jaatharun"] <- par[17])
par[18]/fNS  ## Time when migration between C and H started
(estimates_sp_mod["t start mig C-H", "best.frsp.jaatharun"] <- par[18])
par[19]*fNS*(L-1)  ## recombination rate
(estimates_sp_mod["rho/theta", "best.frsp.jaatharun"] <- par[19]/3.18e-9)
3.18e-9 * fNS * L ## theta (per locus) for the Spanish population

estimates_hd_mod
estimates_sp_mod
