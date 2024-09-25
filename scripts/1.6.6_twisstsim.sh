#!/bin/bash -l
#SBATCH -J twisstsim
#SBATCH --get-user-env
#SBATCH --mail-user=gwee@biologie.uni-muenchen.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=4
#SBATCH --time=2-00:00:00
#SBATCH -o /dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2/slurms/slurm-%j-%x.out

# sbatch 1.6.6_twisstsim.sh 3PopModel8x "-g EURns -g SPA -g EURc -g O --outgroup O" 3popsim run2
# sbatch 1.6.6_twisstsim.sh 3PopModel7x "-g EURns -g SPA -g EURc -g O --outgroup O" 3popsim run2

# sbatch 1.6.6_twisstsim.sh 4PopModel8x "-g cnx3 -g cnx6 -g cor1 -g cor2 -g O --outgroup O" 4popsim run2
# sbatch 1.6.6_twisstsim.sh 4PopModel7x "-g cnx3 -g cnx6 -g cor1 -g cor2 -g O --outgroup O" 4popsim run2

# sbatch 1.6.6_twisstsim.sh 3PopModel8jaatha "-g SPA -g EURc -g EURns -g O --outgroup O" 3popsim run1
# sbatch 1.6.6_twisstsim.sh 3PopModel7jaatha "-g SPA -g EURc -g EURns -g O --outgroup O" 3popsim run1

module load vcftools
module load bcftools
conda activate biotools

ref="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2/genome_HC_allpaths41687_v2.5_chrW.fasta"
dat="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2"
scaff="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/01_probes/scaffold2chr"
gentools="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2/06_results/fst/genomics_general"

echo $(date)
STARTTIME=$(date +%s)

cd ${dat}/fastsimcoal2/4Pop/cor1_cor2to3_cnx1to3_cnx6_50ind_folded_all_nochr18/fastsimcoal2/bestruns/twisst
cd ${1}tree_${4}

python ${dat}/06_results/twisst/twisst.py -t ${1}_${4}.trees.gz -w ${1}_${4}.weights.csv.gz --method complete --groupsFile ${3}.tsv $2

conda activate twisst
Rscript ../plot_twisstsim.R $1 $4
# Rscript ../plot_twisstsim_4pop.R $1 $4
# Rscript ../twisstntern.R $1 $4

ENDTIME=$(date +%s)
echo $(date)
echo "It takes $(($ENDTIME - $STARTTIME)) seconds to complete this task"
echo "It takes $((($ENDTIME - $STARTTIME)/60)) mins to complete this task"
echo "It takes $((($ENDTIME - $STARTTIME)/3600)) hours to complete this task"
