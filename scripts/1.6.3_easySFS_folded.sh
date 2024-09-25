#!/bin/bash -l
#SBATCH -J foldedSFS
#SBATCH --get-user-env
#SBATCH --mail-user=gwee@biologie.uni-muenchen.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=2-00:00:00
#SBATCH -o /dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2/slurms/slurm-%j-%x.out

# sbatch 1.6.3_easySFS_folded.sh 05.1_recal/overlap/neuall_incIRQ/134inds_overlapped_filtered_norepeats_hwe_neuall_incIRQ cor1_cor2to3_cnx1to3_cnx6_50ind_folded_all_nochr18 30,30,30,10 4PopModel1

conda activate easySFS

ref="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2/genome_HC_allpaths41687_v2.5_chrW.fasta"
dat="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2"

echo $(date)
STARTTIME=$(date +%s)

popmodel=${4}
popgrp=${popmodel:0:4}
echo ${popgrp}
poptype=${popmodel:0:1}
echo $poptype

cd ${dat}/fastsimcoal2/${popgrp}

#if folded SFS
python3 ${dat}/fastsimcoal2/easySFS/easySFS.py -i ${dat}/${1}.vcf.gz -p ${2}.pop -a -f --proj ${3} -o $2 --prefix ${4}

cd ./${2}/fastsimcoal2
for i in {2..8}
do
cp 4PopModel1_jointMAFpop1_0.obs 4PopModel${i}_jointMAFpop1_0.obs
cp 4PopModel1_jointMAFpop2_0.obs 4PopModel${i}_jointMAFpop2_0.obs
cp 4PopModel1_jointMAFpop3_0.obs 4PopModel${i}_jointMAFpop3_0.obs
cp 4PopModel1_jointMAFpop2_1.obs 4PopModel${i}_jointMAFpop2_1.obs
cp 4PopModel1_jointMAFpop3_1.obs 4PopModel${i}_jointMAFpop3_1.obs
cp 4PopModel1_jointMAFpop3_2.obs 4PopModel${i}_jointMAFpop3_2.obs
done

ENDTIME=$(date +%s)
echo $(date)
echo "It takes $(($ENDTIME - $STARTTIME)) seconds to complete this task"
echo "It takes $((($ENDTIME - $STARTTIME)/60)) mins to complete this task"
echo "It takes $((($ENDTIME - $STARTTIME)/3600)) hours to complete this task"
