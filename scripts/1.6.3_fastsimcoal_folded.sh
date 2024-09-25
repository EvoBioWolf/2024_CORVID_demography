#!/bin/bash -l
#SBATCH -J 78xfolded
#SBATCH --get-user-env
#SBATCH --mail-user=gwee@biologie.uni-muenchen.de
#SBATCH --clusters=serial
#SBATCH --partition=serial_std
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --time=4-00:00:00
#SBATCH -o /dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2/slurms/slurm-%j-%x.out

# for bs in {1..100};do for i in {7..8}; do sbatch 1.6.3_fastsimcoal_folded.sh cor1_cor2to3_cnx1to3_cnx6_50ind_folded_all_nochr18 30,30,30,10 4PopModel${i}x ${bs}; done; done #serial

conda activate easySFS

ref="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2/genome_HC_allpaths41687_v2.5_chrW.fasta"
dat="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2"

echo $(date)
STARTTIME=$(date +%s)

popmodel=${3}
popgrp=${popmodel:0:4}
echo ${popgrp}
echo run${bs}
cd ${dat}/fastsimcoal2/${popgrp}
cd ${1}/fastsimcoal2

coalsim=1000000 #ideally between 200000 aand 1000000
ECM=100 # at least 20, better between 50 and 100
echo $1 $3

mkdir run${4}
cp ${3}.tpl ${3}.est ${3}_joint*AFpop*.obs run${4}
cd run${4}
${dat}/fastsimcoal2/fsc27 -t ${3}.tpl -e ${3}.est -m -M -L ${ECM} -n ${coalsim} -q -c 12 --foldedSFS


ENDTIME=$(date +%s)
echo $(date)
echo "It takes $(($ENDTIME - $STARTTIME)) seconds to complete this task"
echo "It takes $((($ENDTIME - $STARTTIME)/60)) mins to complete this task"
echo "It takes $((($ENDTIME - $STARTTIME)/3600)) hours to complete this task"
