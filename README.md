# 2024_CORVUS_demography

See [pipeline.md](./pipeline.md) for a detailed documentation of all analyses done, except Jaatha, which is documented in a separate markdown. 
All scripts are stored in the folder [scripts](./scripts/). Intermediate files (such as the input files for [Fastsimcoal](./scripts/fastsimcoal/)) are stored in respective subfolders with the scripts

All analyses are conducted using LRZ BioHPC, LRZ CM2, and LRZ AI clusters. Pipeline written by Chyiyin Gwee (contact chyiyingwee[@]gmail.com)

[pipeline.md](./pipeline.md) in the following order (link to scripts):

* Raw data processing: [trim](./scripts/1.1.1_rawreads_pro.sh), [map ](./scripts/1.2.1_map.sh), [rmdup](./scripts/1.3.1_markdup.sh), [mapDamage](./scripts/1.3.2_mapdamage.sh), [ploidy](./scripts/1.3.4_ploidy.sh)
* Variant calling: [ANGSD](./scripts/1.4.1_angsdrecal.sh), [SAMtools](./scripts/samtools), [GATK](./scripts/1.4.3_GATKrawHap.sh), [overlap](./scripts/1.4.5_overlap.sh), [paralog](./scripts/1.4.5_paralogs.sh), [pop files](./scripts/05.1_recal/) 
* [Combine invariant sites and the final set of variant sites](./scripts/1.4.5_gvcf_combineinvariant.sh)
* [PCA](./scripts/pca.R)
* [Admixture](./scripts/1.5.1_admixture.sh): [plot](./scripts/admixtureplot.R)
* [Admixtools2](./scripts/1.5.3_admixtools.sh): [input files](./scripts/06_results/admixtools/)
* [Dsuite](./scripts/1.5.4_dsuite.sh):  [input files](./scripts/06_results/dsuite/)
* Summary statistics: [TajD](./scripts/1.5.5_summarystats.sh), [FST, Dxy, Da, Pi](./scripts/1.5.6_fst.sh), [plot](./scripts/cmplot.R), [pop files](./scripts/06_results/sumstats/)
* Neutral sites: [variant](./scripts/1.6.1_neutral.sh), [invariant](./scripts/1.6.1_neutral_invariantsites.sh)
* [MSMC2](./scripts/2.1.2_msmc_ind.sh): [input](./scripts/2.1.1_SNPable.sh)
* [Stairway plot](./scripts/2.2.2_stairway_rescaled.sh): [SFS](./scripts/2.2.1_1dsfs.sh), [blueprints](./scripts/stairway/)
* [SMC++](./scripts/2.3.1_smcpp_neuperpopvcf.sh)
* Fastsimcoal2: [SFS](./scripts/1.6.3_easySFS_folded.sh), [parameter estimation](./scripts/1.6.3_fastsimcoal_folded.sh), [plot](./scripts/1.6.4_plot.sh), [simulate SFS from estimates](./scripts/1.6.4_modelfit_folded.sh), [bootstrap](./scripts/1.6.5_simparabs_folded.sh), [input files](./scripts/fastsimcoal/4Pop)
* [TWISST](./scripts/1.5.7_twisst.sh): [input files](./scripts/06_results/twisst/)
* Fastsimcoal to TWISST: [simulate trees](./scripts/1.6.6_fastsimcoalTwisst.sh), [twisst](./scripts/1.6.6_twisstsim.sh), [input files](./scripts/fastsimcoal/4Pop/cor1_cor2to3_cnx1to3_cnx6_50ind_folded_all_nochr18/fastsimcoal2/twisst/)
