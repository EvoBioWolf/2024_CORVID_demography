#!/bin/bash -l
#SBATCH -J fsimplot
#SBATCH --get-user-env
#SBATCH --mail-user=gwee@biologie.uni-muenchen.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=1:00:00
#SBATCH -o /dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2/slurms/slurm-%j-%x.out

# sbatch 1.6.4_plot.sh cor1_cor2to3_cnx1to3_cnx6_50ind_folded_all_nochr18 4PopModel7x 100 cor1 cor2to3 cnx1to3 cnx6
# sbatch 1.6.4_plot.sh cor1_cor2to3_cnx1to3_cnx6_50ind_folded_all_nochr18 4PopModel8x 100 cor1 cor2to3 cnx1to3 cnx6

conda activate easySFS

ref="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2/genome_HC_allpaths41687_v2.5_chrW.fasta"
dat="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2"

echo $(date)
STARTTIME=$(date +%s)

popmodel=${2}
popgrp=${popmodel:0:4}
echo ${popgrp}

cd ${dat}/fastsimcoal2/${popgrp}
cd ${1}/fastsimcoal2

#create dummy files for runs which failed
for i in {1..100}
do
if [ $(ls ./run$i/${2}/$2.bestlhoods | wc -l) -eq 1 ]
then
    echo "file already exist"
else
    echo "generate dummy file" 
    echo "9 0" > ./run$i/${2}/$2.bestlhoods  
fi
done

#smallest difference between observed and expected likelihoods
cat run{1..100}/${2}/${2}.bestlhoods | grep -v "MaxObsLhood\|MaxEstLhood" | awk '{print NR,$(NF-1),$NF}' | awk '{print NR,$3-$2}' | sed -e "s/-9/NA/g" | grep -v "NA" | sort -k2 > ${2}_${3}bs.bestlhoods
cat ./run{1..100}/${2}/${2}.bestlhoods | grep -v NPOP1 | nl | grep -v "9 0" > ${2}_${3}.bestlhoods.all

mkdir bestruns
for i in $(head -1 ${2}_${3}bs.bestlhoods | awk '{print $1}')
do 
cp -r run${i}/${2} ./bestruns
cp run${i}/${2}.est ./bestruns/${2}
cp run${i}/${2}_*.obs run${i}/${2}.tpl run${i}/${2}.par ./bestruns/${2}
done

zless *.bestlhoods > ./bestruns/bestlhoods.txt
cd bestruns
conda activate r4.2
Rscript ${dat}/fastsimcoal2/calculateAIC.sh $2
Rscript ${dat}/fastsimcoal2/SFStools.R -t print2D -i $2 -z
Rscript ${dat}/fastsimcoal2/SFStools.R -t 2Dto1D -i $2 -z
Rscript ${dat}/fastsimcoal2/SFStools.R -t print1D -i $2 -z
Rscript ${dat}/fastsimcoal2/plotModel.R -p $2 -l ${5},${6},${7},${8}

#see which model with lowest AIC
zless ./*/*.AIC > AIC.txt

ENDTIME=$(date +%s)
echo $(date)
echo "It takes $(($ENDTIME - $STARTTIME)) seconds to complete this task"
echo "It takes $((($ENDTIME - $STARTTIME)/60)) mins to complete this task"
echo "It takes $((($ENDTIME - $STARTTIME)/3600)) hours to complete this task"
