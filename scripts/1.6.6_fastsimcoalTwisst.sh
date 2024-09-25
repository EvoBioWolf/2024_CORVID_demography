#!/bin/bash -l
#SBATCH -J twisstfsimcoal
#SBATCH --get-user-env
#SBATCH --mail-user=gwee@biologie.uni-muenchen.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=2-00:00:00
#SBATCH -o /dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2/slurms/slurm-%j-%x.out

# sbatch 1.6.6_fastsimcoalTwisst.sh 3PopModel7x 3popsim run1
# sbatch 1.6.6_fastsimcoalTwisst.sh 3PopModel8x 3popsim run1
# sbatch 1.6.6_fastsimcoalTwisst.sh 4PopModel7x 4popsim run1
# sbatch 1.6.6_fastsimcoalTwisst.sh 4PopModel8x 4popsim run1
# sbatch 1.6.6_fastsimcoalTwisst.sh 3PopModel7x 3popsim run2
# sbatch 1.6.6_fastsimcoalTwisst.sh 3PopModel8x 3popsim run2
# sbatch 1.6.6_fastsimcoalTwisst.sh 4PopModel7x 4popsim run2
# sbatch 1.6.6_fastsimcoalTwisst.sh 4PopModel8x 4popsim run2
# sbatch 1.6.6_fastsimcoalTwisst.sh 3PopModel7jaatha 3popsim run1
# sbatch 1.6.6_fastsimcoalTwisst.sh 3PopModel8jaatha 3popsim run1

module load vcftools/0.1.14-gcc8
module load bcftools/1.10.2-gcc8
conda activate easySFS

ref="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2/genome_HC_allpaths41687_v2.5_chrW.fasta"
dat="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2"

echo $(date)
STARTTIME=$(date +%s)

cd ${dat}/fastsimcoal2/4Pop/cor1_cor2to3_cnx1to3_cnx6_50ind_folded_all_nochr18/fastsimcoal2/bestruns/twisst
#copy the bestrun models and model the parL file for simulation 

#produce treefile
${dat}/fastsimcoal2/fsc27 -i ${1}tree_${3}.par -n 1 -T 2

cd ${1}tree_${3}
cp ../${2}.tsv ./

#modify simulated treefile
awk 'NR>3 { print }' *_true_trees.trees | sed "s/\ttree NumGen_tree_1_1_pos_0 = \[&U\] //g" | gzip > ${1}_${3}.trees.gz

#create windows file
grep "on chromosome" ${1}tree_${3}_1_1.arp | awk '{print $7, $2}' > ${1}_scaffold_sites.win.tmp
grep -A1 "on chromosome" ${1}tree_${3}_1_1.arp | grep -v "on chromosome" | sed "s/#\|,//g" | awk '{print $1, $NF, int(($1+$NF)/2)}' > ${1}_pos.win.tmp
paste ${1}_scaffold_sites.win.tmp ${1}_pos.win.tmp | awk -v OFS="\t" '{print $1, $3, $4, $5, $2, "-100"}' | awk -v OFS="\t" 'BEGIN{print "scaffold", "start", "end", "mid", "sites", "lnL"}1' > ${1}_${3}.win.data.tsv
rm *.tmp

ENDTIME=$(date +%s)
echo $(date)
echo "It takes $(($ENDTIME - $STARTTIME)) seconds to complete this task"
echo "It takes $((($ENDTIME - $STARTTIME)/60)) mins to complete this task"
echo "It takes $((($ENDTIME - $STARTTIME)/3600)) hours to complete this task"
