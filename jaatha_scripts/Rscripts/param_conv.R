## Model 7 and 8, that is with population growth in the French-German population and in the Hooded Crow population

## convert fastsimcoal to msprime frsp parameters

conv_fsc2mspr_sp <- function(p) {

    c(p["NANC3"]/2, p["TDIV3"], p["NANC1"]/2, p["NANC2"]/2, p["TDIV2"], p["NPOP3"]/2,
      -p["GROW2"], p["NPOP4"]/2, p["TDIV1"], p["NPOP2"]/2, p["NPOP1"]/2,
      p["MIG10R"], p["MIG01R"], p["MIG12R"], p["MIG21R"],
      p["MIG32R"], p["MIG23R"], p["TMRMG"], p[19], p["TEGR2"], p["TEGR1"], -p["GROW1"])
}


## convert msprime to fastsimcoal frsp paramters

conv_mspr2fsc_sp <- function(p) {

    c(NANC3=p[1]*2, TDIV3=p[2], NANC1=p[3]*2, NANC2=p[4]*2, TDIV2=p[5],NPOP3=p[6]*2,
      GROW2 = -p[7], NPOP4=p[8]*2, TDIV1=p[9], NPOP2=p[10]*2, NPOP1=p[11]*2,
      MIG10R=p[12], MIG01R=p[13], MIG12R=p[14], MIG21R=p[15],
      MIG32R=p[16], MIG23R=p[17], TMRMG=p[18], p[19], TEGR2=p[20], TEGR1=p[21], GROW1=-p[22])
}

## convert msprime to jaatha data_simulator frsp parameters

conv_mspr2jds_sp <- function(p) {
    
    p[20] <- p[20]/p[5]
    p[21] <- p[21]/p[9]
    p[9] <- p[9] - p[18]
    p[5] <- p[5] - p[9] - p[18]
    p[2]  <- p[2] -  p[5] - p[9] - p[18]
    p[c(1,2,3,4,5,6,8,9,10,11,18,19)] <- log(p[c(1,2,3,4,5,6,8,9,10,11,18,19)])
    p
}

## convert jaatha data_simulator to msprime parameter frsp parameters

conv_jds2mspr_sp <- function(par) {
    ## par is a vector with the following entries:
    ## 1: log of Population size of ancestral population ASI of S and I.
    ## 2: log of Time between the S--I split and the H--I split.
    ## 3: log of Population size of S after this split but before S-W split
    ## 4: log of Population size of I after this split but before H-I split
    ## 5: log of Time between H--I and S--W split
    ## 6: log of Current size if H
    ## 7: Growth rate of H
    ## 8: log of Current size of I (implicitly determins growth rate of I, assuming exponetial change)
    ## 9: log of Time betweeen the split of C from S and the time when migration between C and H started
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

    ## output is a vector with:
    ## 1: Population size of ancestral population ASI of S and I.
    ## 2: Time of the split of S and I.
    ## 3: Population size of S after this split but before S-W split
    ## 4: Population size of I after this split but before H-I split
    ## 5: Time T_H of split between H and I.
    ## 6: Current size if H
    ## 7: Growth rate of H
    ## 8: Current size of I
    ## 9: Time T_W of the split of W from S
    ## 10: Current size of W
    ## 11: Current size of S
    ## 12: Migration rate from S to W (actually trace-back rate of lineages from W back to S)
    ## 13: Migration rate from W to S (actually trace-back rate of lineages from S back to W)
    ## 14: Migration rate from H to W (actually trace-back rate of lineages from W back to H)
    ## 15: Migration rate from W to H (actually trace-back rate of lineages from H back to W)
    ## 16: Migration rate from H to I (actually trace-back rate of lineages from I back to H)
    ## 17: Migration rate from I to H (actually trace-back rate of lineages from H back to I)
    ## 18: Time when migration between W and H started
    ## 19: recombination rate
    ## 20: time when growth of hooded crow population started
    ## 21: time when growth of W population started
    ## 22: growth rate of W
    
    c(exp(par[1]), exp(par[2])+exp(par[5])+exp(par[9])+exp(par[18]), exp(par[3]),
      exp(par[4]), exp(par[5])+exp(par[9])+exp(par[18]), exp(par[6]),
      par[7], exp(par[8]), exp(par[9])+exp(par[18]), exp(par[10]), exp(par[11]),
      par[12], par[13], par[14], par[15], par[16], par[17], exp(par[18]), exp(par[19]),
      par[20]*(exp(par[5])+exp(par[9])+exp(par[18])), par[21]*(exp(par[9])+exp(par[18])),
      par[22])
}

## convert fastsimcoal to msprime frhd parameters

conv_fsc2mspr_hd <- function(p) {

    c(p["NANC3"]/2, p["TDIV3"], p["NANC1"]/2, p["NANC2"]/2, p["TDIV2"], p["NPOP3"]/2,
            -p["GROW2"], p["NPOP4"]/2, p["TDIV1"], p["NPOP2"]/2, p["NPOP1"]/2,
            p["MIG10R"], p["MIG01R"], p["MIG12R"], p["MIG21R"],
            p["MIG32R"], p["MIG23R"], p["TMRMG"], p[19], p["TEGR2"], p["TEGR1"], -p["GROW1"])
}

## convert msprime to fastsimcoal frhd paramters

conv_mspr2fsc_hd <- function(p) {

    c(NANC3=p[1]*2, TDIV3=p[2], NANC1=p[3]*2, NANC2=p[4]*2, TDIV2=p[5],NPOP3=p[6]*2,
            GROW2 = -p[7], NPOP4=p[8]*2, TDIV1=p[9], NPOP2=p[10]*2, NPOP1=p[11]*2,
            MIG10R=p[12], MIG01R=p[13], MIG12R=p[14], MIG21R=p[15],
            MIG32R=p[16], MIG23R=p[17], TMRMG=p[18], p[19],  TEGR2=p[20], TEGR1=p[21], GROW1=-p[22])
}

## convert msprime to jaatha data_simulator frhd parameters

conv_mspr2jds_hd <- function(p) {

    p[20] <- p[20]/p[9]
    p[21] <- p[21]/p[9]
    p[9] <- p[9] - p[18]
    p[5] <- p[5] - p[9] - p[18]
    p[2]  <- p[2] - p[5] - p[9] - p[18]
    p[c(1,2,3,4,5,6,8,9,10,11,18,19)]  <- log(p[c(1,2,3,4,5,6,8,9,10,11,18,19)])
    p
}

## convert jaatha data_simulator to msprime parameter frhd parameters

conv_jds2mspr_hd <- function(par) {
   
    ## entries of input par are:
    ## 1: log of Population size of ancestral population SI of S and I.
    ## 2: log of Time between the S--I split and the H--I split.
    ## 3: log of size of HC at the time of the H--C split
    ## 4: log of Population size of HI, that is I after this split but before H-I split
    ## 5: log of Time between H--I and H--C split
    ## 6: log of Population size of H at present time
    ## 7: Growth rate of H (and HI) since H--I split
    ## 8: log of Current size of I
    ## 9: log of Time betweeen the split of C from H and the time when migration between C and S started
    ## 10: log of Current size of C
    ## 11: log of Current size of S
    ## 12: Migration rate from S to C (actually trace-back rate of lineages from C back to S)
    ## 13: Migration rate from C to S (actually trace-back rate of lineages from S back to C)
    ## 14: Migration rate from H to C (actually trace-back rate of lineages from C back to H)
    ## 15: Migration rate from C to H (actually trace-back rate of lineages from H back to C)
    ## 16: Migration rate from H to I (actually trace-back rate of lineages from I back to H)
    ## 17: Migration rate from I to H (actually trace-back rate of lineages from H back to I)
    ## 18: log of the time when the migration between C and S started
    ## 19: log of recombination rate
    ## 20: fraction of time of growth of hooded crow population compared to existence since split from W
    ## 21: fraction of time of growth of W population relative to existence since split from S
    ## 22: growth rate of W

    ## output:
    ## 1: Population size of ancestral population ASI of S and I.
    ## 2: Time of the split of S and I.
    ## 3: Population size of H right before H-W split
    ## 4: Population size of H-I ancestral after this split but before H-I split
    ## 5: Time T_H of split between H and I.
    ## 6: Current size if H
    ## 7: Growth rate of H since split (implicitly determines size of H at time T_H)
    ## 8: Current size of I
    ## 9: Time of the split of W from H
    ## 10: Current size of W
    ## 11: Current size of S
    ## 12: Migration rate from S to W (actually trace-back rate of lineages from W back to S)
    ## 13: Migration rate from W to S (actually trace-back rate of lineages from S back to W)
    ## 14: Migration rate from H to W (actually trace-back rate of lineages from W back to H)
    ## 15: Migration rate from W to H (actually trace-back rate of lineages from H back to W)
    ## 16: Migration rate from H to I (actually trace-back rate of lineages from I back to H)
    ## 17: Migration rate from I to H (actually trace-back rate of lineages from H back to I)
    ## 18: Time when migration between W and S started
    ## 19: recombination rate
    ## 20: time when growth of hooded crow population started
    ## 21: time when growth of W population started
    ## 22: growth rate of W
    
    c(exp(par[1]), exp(par[2])+exp(par[5])+exp(par[9])+exp(par[18]), exp(par[3]), exp(par[4]),
      exp(par[5])+exp(par[9])+exp(par[18]), exp(par[6]),
      par[7], exp(par[8]), exp(par[9])+exp(par[18]), exp(par[10]), exp(par[11]),
      par[12], par[13], par[14], par[15], par[16], par[17], exp(par[18]), exp(par[19]),
      par[20]*(exp(par[9])+exp(par[18])), par[21]*(exp(par[9])+exp(par[18])),
      par[22])
}

if(FALSE) {
    
    all(round(conv_mspr2jds_sp(conv_jds2mspr_sp(1:22)),8) == 1:22)
    all(round(conv_mspr2jds_hd(conv_jds2mspr_hd(1:22)),8) == 1:22)
    all(round(conv_fsc2mspr_hd(conv_mspr2fsc_hd(1:22)),8) == 1:22)
    all(round(conv_fsc2mspr_sp(conv_mspr2fsc_sp(1:22)),8) == 1:22)
}
