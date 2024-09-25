# 2024_CORVID_demography

See `pipeline.md` for all analyses done, except Jaatha, which is documented in a separate markdown. 
All scripts are stored in the folder `scripts`. Intermediate files (such as the `tpl` and `est` input files for Fastsimcoal) are stored in respective subfolders with the scripts

All analyses are conducted using LRZ BioHPC, LRZ CM2, and LRZ AI clusters. Pipeline written by Chyiyin Gwee (contact chyiyingwee[@]gmail.com)

`pipeline.md` in the following order (with the corresponding scripts in bracket):

`tip: use sidebar to skip to a particular analysis of the markdown` 

* Raw data processing (1.1*-1.3*)
* Variant calling (1.4*)
* Combine invariant sites and the final set of variant sites (1.4.5*)
* PCA (pca.R)
* Admixture (1.5.1*, admixtureplot.R)
* Admixtools2 (1.5.3*)
* Dsuite (1.5.4*)
* Summary statistics: FST, Dxy, Da, Pi (1.5.5*-1.5.6*, cmplot.R)
* Neutral sites (1.6.1*)
* MSMC2 (2.1*)
* Stairway plot (2.2*)
* SMC++ (2.3*)
* Fastsimcoal2 (1.6.3*-1.6.5*)
* TWISST (1.5.7*)
* Fastsimcoal to TWISST (1.6.6*)
