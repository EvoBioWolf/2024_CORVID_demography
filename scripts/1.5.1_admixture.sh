#!/bin/bash -l
#SBATCH -J admix
#SBATCH --get-user-env
#SBATCH --mail-user=gwee@biologie.uni-muenchen.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_production
#SBATCH --cpus-per-task=4
#SBATCH --time=7-00:00:00
#SBATCH -o /dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2/slurms/slurm-%j-%x.out
#SBATCH --mem-per-cpu=4763mb

# for i in {1..10}; do sbatch 1.5.1_admixture.sh 134inds_overlapped_filtered_norepeats_ldpruned $i poplist_134 overlap; done

conda activate biotools
#picard.jar=2.25.7

ref="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2/genome_HC_allpaths41687_v2.5_chrW.fasta"
dat="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2"

echo $(date)
STARTTIME=$(date +%s)

cd $dat/05.1_recal/$4
admixture --cv ${1}.bed $2 -j4 | tee log${2}.out 


cd ${dat}/06_results
mkdir admixture
cd admixture
mkdir $1
cd ./$1
mv $dat/05.1_recal/$4/log${2}.out ./
mv $dat/05.1_recal/$4/${1}.${2}.Q ./
mv $dat/05.1_recal/$4/${1}.${2}.P ./

conda activate renv
Rscript ${dat}/admixtureplot.R $1 $2 ${dat}/06_results/admixture/$1 ../$3.txt

ENDTIME=$(date +%s)
echo $(date)
echo "It takes $(($ENDTIME - $STARTTIME)) seconds to complete this task"
echo "It takes $((($ENDTIME - $STARTTIME)/60)) mins to complete this task"
echo "It takes $((($ENDTIME - $STARTTIME)/3600)) hours to complete this task"


