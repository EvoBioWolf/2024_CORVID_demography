#!/bin/bash -l
#SBATCH -J dsuite
#SBATCH --get-user-env
#SBATCH --mail-user=gwee@biologie.uni-muenchen.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=4
#SBATCH --time=2-00:00:00
#SBATCH -o /dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2/slurms/slurm-%j-%x.out
#SBATCH --mem-per-cpu=4763mb

# sbatch 1.5.4_dsuite.sh 05.1_recal/overlap/134inds_overlapped_filtered_norepeats_hwe_AMcrow_biallele fbranch

conda activate biotools
#picard.jar=2.25.7

ref="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2/genome_HC_allpaths41687_v2.5_chrW.fasta"
dat="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2"
vcf="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2/05.1_recal/overlap/134inds_overlapped_filtered_norepeats_outgroup.vcf.gz"
scaff="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/01_probes/scaffold2chr"
out="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2/06_results/dsuite/dinvestigate"

echo $(date)
STARTTIME=$(date +%s)

cd $dat/06_results
mkdir dsuite
cd dsuite

./Dsuite/Build/Dsuite Dtrios ${dat}/${1}.vcf.gz SETS_$2.txt -t tree.nwk
./Dsuite/Build/Dsuite Dinvestigate ${dat}/${1}.vcf.gz SETS_$2.txt SETS_trios2.txt -n $2
mv *_$2_50_25.txt ./dinvestigate

./Dsuite/Build/Dsuite Fbranch tree.nwk SETS_$2_tree.txt > fbranch_matrix.txt
#plotting require python3.8
conda activate py3.8
python ./Dsuite/utils/dtools.py fbranch_matrix.txt tree.nwk

# cd dinvestigate
#mv cnx6_cnx*_cor*_50_25.txt ./cnx6_cnx_cor
#mv cnx*_cnx6_pec1_50_25.txt ./cnx_cnx6_pec1
#mv cor*_cor*_cnx*_50_25.txt ./cor_cor_cnx
#mv ori*_ori*_cnx*_50_25.txt ./ori_ori_cnx
#mv ori*_ori*_pec1_50_25.txt ./ori_ori_pec1


ENDTIME=$(date +%s)
echo $(date)
echo "It takes $(($ENDTIME - $STARTTIME)) seconds to complete this task"
echo "It takes $((($ENDTIME - $STARTTIME)/60)) mins to complete this task"
echo "It takes $((($ENDTIME - $STARTTIME)/3600)) hours to complete this task"


