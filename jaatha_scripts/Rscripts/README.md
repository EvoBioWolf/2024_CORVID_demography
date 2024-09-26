The R scripts `jaatha_frsp_scrm_folded.R` and `jaatha_frhd_scrm_folded.R`
apply jaatha to estimate the parameters of the demographic models in which
the Western European population first split off from the south-western or
from the hooded crow population, respectively.

For further analyses as documented below, the jaatha output `jfrsp_zi_123_folded.RData`
should be moved to a folder `jaatha_estimations`.

`eval_parbs.R` contains R code to explore parameter estimations and compare likelihoods.

The R scripts `finite_sites_simul.R` and `is_simul.R` can be used to simulate data
according to the two demographic models with given parameter values.
While `finite_sites_simul.R` assumes a finite-sites mutation model, `is_simul.R`
simulates data based on infinite-sites assumptions.
For subsequent analyses the scripts documented below assume that the simulated data
should are moved to folders `fs_sim_folded` and `is_sim_folded`, respectively.

`compare_jsfs.R` can be used to compare the pairwise joint frequency spectra of
simulated data to those of the empirical data.

`jaatha_frsp_parbs_folded.R` carries out parametric-bootstrap re-estimations
using simulted paramtric-bootstrap simulations in folder `fs_sim_folded`.
The re-estimations can then be moved to a folder `parbs_frsp_folded` and be
evaluated with the R script `eval_parbs.R`.

The R scripts `sumstat_funcs.R`, `param_conv.R` and `jaatha_frsp_scrm_func.R` are
used by the other R files and contain functions to calculate the summary statistics
used in the jaatha analyses, to convert parameters according to the different scalings
used by `fastsimcoal`, `msprime` and `jaatha` and to simulate coalescent trees with
recombination according to the demographic models.
