#!/bin/bash -l
#SBATCH -J admixall
#SBATCH --get-user-env
#SBATCH --mail-user=gwee@biologie.uni-muenchen.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_production
#SBATCH --cpus-per-task=4
#SBATCH --ntasks=1
#SBATCH --time=7-00:00:00
#SBATCH -o /dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2/slurms/slurm-%j-%x.out

# sbatch 1.5.3_admixtools.sh 134inds_overlapped_filtered_norepeats_hwe_outgroup_biallele 134inds_all poplist_134 1 100 2 01_134inds_outgroup_all
# sbatch 1.5.3_admixtools.sh 134inds_overlapped_filtered_norepeats_hwe_AMcrow_biallele AM_134inds_all poplist_134_AM 1 100 2 01_134inds_AM_all

module load vcftools

#dat="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/02_outgroups/08_corvuscombine"
#pop="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/demo/poplist"

dat="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2"

echo $(date)
STARTTIME=$(date +%s)

cd $dat/06_results
mkdir admixtools
cd admixtools

#make eigenstrat input
python vcf2eigenstrat.py -v ${dat}/05.1_recal/overlap/${1}.vcf.gz -o $2

rm ${2}.ind
cp ${3}.ind ${2}.ind
awk '{print $2":"$4, $2, $3, $4, $5, $6}' ${2}.snp | awk '{gsub("scaffold_", "", $2)}1' | awk '{print $1, $2+1, $4*0.000001, $4, $5, $6}' OFS="\t" > ${2}_final.snp
mv ${2}_final.snp ${2}.snp

conda activate r4.2
cd ${dat}/06_results/admixtools
Rscript ${dat}/06_results/admixtools/admixtools_1000gen.R $2 $6 $4 $5

#summary & compile
mkdir ${7}
mv ${2}*_adm${6}* ./${7}
cd ${7} 
for i in *_bestscores.txt; do echo $i >> 01_summary.txt; awk 'NR>1 {print $2, $1}' $i | sort -n | head -1 >> 01_summary.txt; done
mkdir 01_bestgraphs
for i in *_bestscores.txt; do base=${i%_bs1to100_bestscores.txt*}; sort -k2 -n $i | awk 'NR>1 {print $1}' | cp $(echo ${base}_run_$(head -n 1).tiff) ./01_bestgraphs; done
conda activate base
cd 01_bestgraphs
for i in *.tiff; do base=${i%_adm*}; convert -append ${base}_*.tiff merged_${base}.tiff; done

ENDTIME=$(date +%s)
echo $(date)
echo "It takes $(($ENDTIME - $STARTTIME)) seconds to complete this task"
echo "It takes $((($ENDTIME - $STARTTIME)/60)) mins to complete this task"
echo "It takes $((($ENDTIME - $STARTTIME)/3600)) hours to complete this task"
