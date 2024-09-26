#!/usr/bin/env python3

## run with python 3.7 or newer!

import numpy as np
import msprime
import sys

from multiallelic import *

NSI = float(sys.argv[1])   ## Population size of ancestral population ASI of S and I.
TSI = float(sys.argv[2])   ## Time of the split of S and I.
NSW = float(sys.argv[3])   ## Population size of S after this split but before S-W split
NHI = float(sys.argv[4])   ## Population size of I after this split but before H-I split
THI = float(sys.argv[5])   ## Time T_H of split between H and I.
NH  = float(sys.argv[6])   ## Current size if H
gH  = float(sys.argv[7])   ## Growth rate of H since split (implicitly determines size of H at time T_H)
NI  = float(sys.argv[8])   ## Current size of I (implicitly determins growth rate of I, assuming exponetial change)
TWS = float(sys.argv[9])   ## Time T_W of the split of W from S
NW  = float(sys.argv[10])  ## Current size of W
NS  = float(sys.argv[11])  ## Current size of S
mWS = float(sys.argv[12])  ## Migration rate from S to W (actually trace-back rate of lineages from W back to S)
mSW = float(sys.argv[13])  ## Migration rate from W to S (actually trace-back rate of lineages from S back to W)
mWH = float(sys.argv[14])  ## Migration rate from H to W (actually trace-back rate of lineages from W back to H)
mHW = float(sys.argv[15])  ## Migration rate from W to H (actually trace-back rate of lineages from H back to W)
mIH = float(sys.argv[16])  ## Migration rate from H to I (actually trace-back rate of lineages from I back to H)
mHI = float(sys.argv[17])  ## Migration rate from I to H (actually trace-back rate of lineages from H back to I)
THW = float(sys.argv[18])  ## Time when migration between W and H started
recomb_rate   = float(sys.argv[19])  ## recombination rate
TgH = float(sys.argv[20])        ## time when growth of hooded crow population started
TgW = float(sys.argv[21])        ## time when growth of W population started
gW  = float(sys.argv[22])        ## growth rate of W
random_seed   = int(float(sys.argv[23]))
folded = sys.argv[24]=="TRUE"

## if values are needed for testing purposes set to True:
testing = False
if testing :
    NSI = 100000
    TSI = 100000 
    NSW = 10000
    NHI = 10000 
    THI = 8000
    NH  = 10000 
    gH  = 0.00001
    NI  = 10000 
    TWS = 1000
    NW  = 10000 
    NS  = 10000
    mWS =  0.001
    mSW =  0.001
    mWH =  0.001
    mHW =  0.001
    mIH =  0.001
    mHI =  0.001
    recomb_rate = 0.00001
    TgH = 2000
    TgW = 800
    gW  = 0.00001
    random_seed = 1239

seqlen = 1_000
nloci2  = round((625134484+11146221)/1000)
## nloci2  = round(1.2e9*11266153/14.823e6/1000)  ## simulations without recombination for frequency spectra
## nloci2  = round(1.2e9*9359926/14.823e6/1000)

demography = msprime.Demography()
demography.add_population(name="S", initial_size = NS)
demography.add_population(name="W", initial_size = NW, growth_rate = gW)
demography.add_population(name="H", initial_size = NH, growth_rate = gH)
demography.add_population(name="I", initial_size = NI)
demography.add_population(name="SW", initial_size = NSW)
demography.add_population(name="HI", initial_size = NHI)
demography.add_population(name="SI", initial_size = NSI)
demography.add_population_parameters_change(time=TgH, growth_rate = 0, population="H")
demography.add_population_parameters_change(time=TgW, growth_rate = 0, population="W")
demography.set_migration_rate(source = "S", dest = "W", rate = mSW)
demography.set_migration_rate(source = "W", dest = "S", rate = mWS)
demography.set_migration_rate(source = "W", dest = "H", rate = mWH)
demography.set_migration_rate(source = "H", dest = "W", rate = mHW)
demography.set_migration_rate(source = "I", dest = "H", rate = mIH)
demography.set_migration_rate(source = "H", dest = "I", rate = mHI)
demography.add_migration_rate_change(time = TWS, rate=0, source = "S", dest = "W")
demography.add_migration_rate_change(time = TWS, rate=0, source = "W", dest = "S")
demography.add_migration_rate_change(time = THW, rate=0, source = "W", dest = "H")
demography.add_migration_rate_change(time = THW, rate=0, source = "H", dest = "W")
demography.add_migration_rate_change(time = THI, rate=0, source = "I", dest = "H")
demography.add_migration_rate_change(time = THI, rate=0, source = "H", dest = "I")
demography.add_population_split(time = TWS, derived=["W"], ancestral="SW")
demography.add_population_split(time = TWS, derived=["S"], ancestral="SW")
demography.add_population_parameters_change(time=TWS, initial_size=NSW, population="SW")
demography.add_population_split(time = THI, derived=["H"], ancestral="HI")
demography.add_population_split(time = THI, derived=["I"], ancestral="HI")
demography.add_population_parameters_change(time=THI, initial_size=NHI, population="HI")
demography.add_population_split(time = TSI, derived=["HI"], ancestral="SI")
demography.add_population_split(time = TSI, derived=["SW"], ancestral="SI")
demography.add_population_parameters_change(time=TSI, initial_size=NSI, population="SI")

demography.sort_events()
## print(demography.debug())

sampsize = (15, 15, 15, 5)  ## number of sampled diploids from individuals from S, W, H and I
    
tss2 = msprime.sim_ancestry(samples = [msprime.SampleSet(sampsize[0], "S"),
                                       msprime.SampleSet(sampsize[1], "W"),
                                       msprime.SampleSet(sampsize[2], "H"),
                                       msprime.SampleSet(sampsize[3], "I")],
                            demography = demography,
                            recombination_rate = 0.0,
                            model=msprime.SmcPrimeApproxCoalescent(),
                            sequence_length = seqlen,
                            random_seed=random_seed+1,
                            num_replicates = nloci2
                            )



## jsfs = {}
## for i in range(3) :
##     for j in range(i+1, 4) :
##        jsfs[(i,j)] = np.zeros((sampsize[i]*2+1, sampsize[j]*2+1))

fjsfs = {}
for i in range(3) :
    for j in range(i+1, 4) :
        fjsfs[(i,j)] = np.zeros((sampsize[i]*2+1, sampsize[j]*2+1))
   

d = 0
for ts in tss2 :
    d += 1
    ## if d % 10_000 == 0:
    ##      print(d)
    mts = msprime.sim_mutations(ts, rate=3.18e-9, model=msprime.HKY(kappa=2))
    if mts.get_num_sites() > 0:
        for i in range(3) :
            for j in range(i+1, 4) :
                ## jsfs[(i,j)] +=  mts.allele_frequency_spectrum(sample_sets = [ts.samples(population = i),
                ##                                                             ts.samples(population = j)],
                ##                      span_normalise = False, polarised = True)
                fjsfs[(i,j)] += filteredjsfs(mts, samplesets = [ts.samples(population = i),
                                                    ts.samples(population = j)])


for i in range(3) :
    for j in range(i+1, 4) :
        print(f'\njoint allele frequency spectrum for populations {demography.populations[i].name} and {demography.populations[j].name}:')
        n, m = fjsfs[(i,j)].shape
        for k in range(n) :
            for l in range(m) :
                print(int(fjsfs[(i,j)][k,l]), end=" ")
            print("")
