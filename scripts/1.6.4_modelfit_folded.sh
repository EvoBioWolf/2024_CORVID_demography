#!/bin/bash -l
#SBATCH -J modelfit
#SBATCH --get-user-env
#SBATCH --mail-user=gwee@biologie.uni-muenchen.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=2-00:00:00
#SBATCH -o /dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2/slurms/slurm-%j-%x.out

# for i in {7..8}; do sbatch 1.6.4_modelfit_folded.sh 4PopModel${i}_fastsimcoal 4PopModel${i} fastsimcoal cor1_cor2to3_cnx1to3_cnx6_50ind_folded_all_nochr18; done
# for i in {7..7}; do sbatch 1.6.4_modelfit_folded.sh 4PopModel${i}_jaatha 4PopModel${i} jaatha cor1_cor2to3_cnx1to3_cnx6_50ind_folded_all_nochr18; done
# for i in {1..2}; do sbatch 1.6.4_modelfit_folded.sh 4PopModel${i}_fastsimcoal 4PopModel${i} fastsimcoal cor1_cor2to3_cnx1to3_cnx6_50ind_folded_all_nochr18; done
# for i in {7..8}; do sbatch 1.6.4_modelfit_folded.sh 4PopModel${i}x_fastsimcoal 4PopModel${i}x fastsimcoal cor1_cor2to3_cnx1to3_cnx6_50ind_folded_all_nochr18; done
# sbatch 1.6.4_modelfit_folded.sh 4PopModel7x_fastsimcoal 4PopModel7x fastsimcoal cor1_cor2to3_cnx1to3_cnx6_50ind_folded_all_nochr18

conda activate easySFS

ref="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2/genome_HC_allpaths41687_v2.5_chrW.fasta"
dat="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2"

echo $(date)
STARTTIME=$(date +%s)

cd ${dat}/fastsimcoal2/4Pop/${4}/fastsimcoal2
cd bestruns/modelfit

#infinite sites
${dat}/fastsimcoal2/fsc27 -i ${1}_infinite.par -j -m -s0 -x -I -q -n100 -c 4 --foldedSFS

#finite sites
${dat}/fastsimcoal2/fsc27 -i ${1}_finite.par -j -m -s0 -x -q -n100 -c 4 --foldedSFS

##legend##
#-j: output one simulated or bootstrapped SFS per file in a separate directory for easier analysis
#-m: Computes the site frequency spectrum (SFS) for the minor alleles for each population sample and for all pairs of samples (joint SFS)
#-s0: output all SNPs in the DNA sequence
#-x: Does not generate Arlequin output file
#-I:Generates DNA mutations according to an infinite site (IS) mutation model
#-q: quiet mode
#-n100: 100 simulations

#remove the removed sites message on finite SFS
cd $1_finite
for i in {1..100}
do
cd $1_finite_${i}
sed -i 's/(.*)//g' $1_finite_jointMAFpop1_0.obs
sed -i 's/(.*)//g' $1_finite_jointMAFpop2_0.obs
sed -i 's/(.*)//g' $1_finite_jointMAFpop3_0.obs
sed -i 's/(.*)//g' $1_finite_jointMAFpop2_1.obs
sed -i 's/(.*)//g' $1_finite_jointMAFpop3_1.obs
sed -i 's/(.*)//g' $1_finite_jointMAFpop3_2.obs
cd ../
done

cd ${dat}/fastsimcoal2/4Pop/${4}/fastsimcoal2/bestruns/modelfit
conda activate r4.2

for site in finite infinite
do 
cd $1_${site}
awk 'NR>1' $1_${site}.lhoodObs | tr '\t' '\n' | awk 'NR>2' | nl | sort -n -k2 -r > $1_${site}.lhoodObs_sorted
for i in $(head -1 $1_${site}.lhoodObs_sorted | awk '{print $1}')
do
mkdir ${2}
cp *_${i}/*obs ./${2}
cd ./${2}
rename _${3}_${site}_ _ *.obs
rename .obs .txt *.obs
cd ../
cp ../${2}*.obs ./${2}
Rscript ${dat}/fastsimcoal2/SFStools.R -t print2D -i $2 -z
done
cd ../
done

ENDTIME=$(date +%s)
echo $(date)
echo "It takes $(($ENDTIME - $STARTTIME)) seconds to complete this task"
echo "It takes $((($ENDTIME - $STARTTIME)/60)) mins to complete this task"
echo "It takes $((($ENDTIME - $STARTTIME)/3600)) hours to complete this task"
