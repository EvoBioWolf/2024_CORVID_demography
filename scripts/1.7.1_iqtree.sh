#!/bin/bash -l
#SBATCH -J iqtree
#SBATCH --get-user-env
#SBATCH --mail-user=gwee@biologie.uni-muenchen.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_production
#SBATCH --cpus-per-task=10
#SBATCH --time=14-00:00:00
#SBATCH -o /dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2/slurms/slurm-%j-%x.out

# sbatch 1.7.1_iqtree.sh 134inds_overlapped_filtered_norepeats
# sbatch 1.7.1_iqtree.sh 134inds_overlapped_filtered_norepeats_hwe  
# sbatch 1.7.1_iqtree.sh 134inds_overlapped_filtered_norepeats_hwe_monedula

module load vcftools
module load bcftools
conda activate biotools
#picard.jar=2.25.7
#module load plink2
#using plink1.9 instead

ref="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2/genome_HC_allpaths41687_v2.5_chrW.fasta"
dat="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2"
scaff="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/01_probes/scaffold2chr"
gentools="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2/06_results/fst/genomics_general"

echo $(date)
STARTTIME=$(date +%s)

cd $dat/06_results/iqtree
python vcf2phylip/vcf2phylip.py -i ${dat}/05.1_recal/overlap/$1.vcf.gz
 
iqtree -s ${1}.min4.phy -m GTR+G -B 1000 -T 10

#-m MFP+ASC cant be used due to invariant sites
#run splittree later
#-m MFP required >134 CPU ram
ENDTIME=$(date +%s)
echo $(date)
echo "It takes $(($ENDTIME - $STARTTIME)) seconds to complete this task"
echo "It takes $((($ENDTIME - $STARTTIME)/60)) mins to complete this task"
echo "It takes $((($ENDTIME - $STARTTIME)/3600)) hours to complete this task"
