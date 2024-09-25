#!/bin/bash -l
#SBATCH -J fst
#SBATCH --get-user-env
#SBATCH --mail-user=gwee@biologie.uni-muenchen.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=8
#SBATCH --time=2-00:00:00
#SBATCH -o /dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2/slurms/slurm-%j-%x.out


# sbatch 1.5.6_fst.sh 134inds_overlapped_filtered_norepeats_hwe_incinvariants cnx3P cor2 hz1 poplist_134
# sbatch 1.5.6_fst.sh 134inds_overlapped_filtered_norepeats_hwe_incinvariants cnx3P cor1 hz1 poplist_134
# sbatch 1.5.6_fst.sh 134inds_overlapped_filtered_norepeats_hwe_incinvariants cnx3P cor3 hz1 poplist_134
# sbatch 1.5.6_fst.sh 134inds_overlapped_filtered_norepeats_hwe_incinvariants cnx6 cor2 hz1 poplist_134
# sbatch 1.5.6_fst.sh 134inds_overlapped_filtered_norepeats_hwe_incinvariants cnx6 cor1 hz1 poplist_134
# sbatch 1.5.6_fst.sh 134inds_overlapped_filtered_norepeats_hwe_incinvariants cnx6 cor3 hz1 poplist_134
# sbatch 1.5.6_fst.sh 134inds_overlapped_filtered_norepeats_hwe_incinvariants cnx1 cor2 hz1 poplist_134
# sbatch 1.5.6_fst.sh 134inds_overlapped_filtered_norepeats_hwe_incinvariants cnx1 cor1 hz1 poplist_134
# sbatch 1.5.6_fst.sh 134inds_overlapped_filtered_norepeats_hwe_incinvariants cnx1 cor3 hz1 poplist_134
# sbatch 1.5.6_fst.sh 134inds_overlapped_filtered_norepeats_hwe_incinvariants cnx2 cor2 hz1 poplist_134
# sbatch 1.5.6_fst.sh 134inds_overlapped_filtered_norepeats_hwe_incinvariants cnx2 cor1 hz1 poplist_134
# sbatch 1.5.6_fst.sh 134inds_overlapped_filtered_norepeats_hwe_incinvariants cnx2 cor3 hz1 poplist_134
# sbatch 1.5.6_fst.sh 134inds_overlapped_filtered_norepeats_hwe_incinvariants cnx3P cnx6 hz1 poplist_134
# sbatch 1.5.6_fst.sh 134inds_overlapped_filtered_norepeats_hwe_incinvariants cnx1 cnx6 hz1 poplist_134
# sbatch 1.5.6_fst.sh 134inds_overlapped_filtered_norepeats_hwe_incinvariants cnx2 cnx6 hz1 poplist_134

# sbatch 1.5.6_fst.sh 134inds_overlapped_filtered_norepeats_hwe_incinvariants cnx3S cor1 hz1 poplist_134
# sbatch 1.5.6_fst.sh 134inds_overlapped_filtered_norepeats_hwe_incinvariants cnx3S cor2 hz1 poplist_134
# sbatch 1.5.6_fst.sh 134inds_overlapped_filtered_norepeats_hwe_incinvariants cnx3S cor3 hz1 poplist_134
# sbatch 1.5.6_fst.sh 134inds_overlapped_filtered_norepeats_hwe_incinvariants cnx4 cor1 hz1 poplist_134
# sbatch 1.5.6_fst.sh 134inds_overlapped_filtered_norepeats_hwe_incinvariants cnx4 cor2 hz1 poplist_134
# sbatch 1.5.6_fst.sh 134inds_overlapped_filtered_norepeats_hwe_incinvariants cnx4 cor3 hz1 poplist_134
# sbatch 1.5.6_fst.sh 134inds_overlapped_filtered_norepeats_hwe_incinvariants cnx5 cor1 hz1 poplist_134
# sbatch 1.5.6_fst.sh 134inds_overlapped_filtered_norepeats_hwe_incinvariants cnx5 cor2 hz1 poplist_134
# sbatch 1.5.6_fst.sh 134inds_overlapped_filtered_norepeats_hwe_incinvariants cnx5 cor3 hz1 poplist_134

# sbatch 1.5.6_fst.sh 134inds_overlapped_filtered_norepeats_hwe_incinvariants ori1 cnx6 hz1 poplist_134
# sbatch 1.5.6_fst.sh 134inds_overlapped_filtered_norepeats_hwe_incinvariants ori2 cnx6 hz1 poplist_134
# sbatch 1.5.6_fst.sh 134inds_overlapped_filtered_norepeats_hwe_incinvariants ori3 cnx6 hz1 poplist_134
# sbatch 1.5.6_fst.sh 134inds_overlapped_filtered_norepeats_hwe_incinvariants pec1 cnx6 hz1 poplist_134


module load vcftools
module load bcftools
conda activate biotools
#picard.jar=2.25.7
#module load plink2
#using plink1.9 instead

ref="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2/genome_HC_allpaths41687_v2.5_chrW.fasta"
dat="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2"
scaff="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/01_probes/scaffold2chr"

echo $(date)
STARTTIME=$(date +%s)

cd $dat/06_results
mkdir fst
cd fst

if [ $(ls ${1}.geno.gz | wc -l) -eq 1 ]
then
    echo "file already exist, skip this step"
else
    echo "generate geno file" 
    python genomics_general/VCF_processing/parseVCF.py -i ${dat}/05.1_recal/overlap/${1}.vcf.gz --skipIndels --ploidyMismatchToMissing --ploidy 2 -o ${1}.geno.gz
fi

#Simon Martin's script 
python genomics_general/popgenWindows.py --popsFile ${5}.txt -g ${1}.geno.gz \
-o ${4}_${2}_${3}.win.csv \
-f phased -w 50000 -p $2 -p $3 -T 4 --ploidy 2

awk -F ',' '{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $1, $1":"$4}' ${4}_${2}_${3}.win.csv | sed '0,/scaffold/s//CHROM/' | nl > ${4}_${2}_${3}.win.chr.tmp
#1  CHROM start end mid sites pi_$2 pi_$3 dxy_$2_$3 Fst_$2_$3 scaffold scaffold:mid

cd $scaff
for fname in *.scaffolds
do
base=${fname%.scaffolds*}
for i in $(cat $base.scaffolds)
do
echo $i $base
awk -v scaff=$i -v chr=$base '{gsub("^"scaff"$", chr, $2)} 1' ${dat}/06_results/fst/${4}_${2}_${3}.win.chr.tmp > ${dat}/06_results/fst/${4}_${2}_${3}.tmp && mv \
-v ${dat}/06_results/fst/${4}_${2}_${3}.tmp ${dat}/06_results/fst/${4}_${2}_${3}.win.chr.tmp
done
done

#calculate dA
cd $dat/06_results/fst
awk 'NR>1 {print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $9-($7+$8)/2, $11, $12}' ${4}_${2}_${3}.win.chr.tmp | awk 'BEGIN{print "1", "CHROM", "start", "end", "mid", "sites", "pi_pop1", "pi_pop2", "dxy", "Fst", "da", "scaffold", "scaffold:mid"}1' > ${4}_${2}_${3}.win.chr
grep chr18 ${4}_${2}_${3}.win.chr | awk 'BEGIN{print "1", "CHROM", "start", "end", "mid", "sites", "pi_pop1", "pi_pop2", "dxy", "Fst", "da", "scaffold", "scaffold:mid"}1' > ${4}_${2}_${3}.win.chr18
conda activate renv
Rscript ${dat}/cmplot.R ${2}_${3} ${4}

#non-elevated chr18 sites
for i in {1..3}
do
grep -v $(awk -v ORS="\\\|" 1  elevatedregions.pos | awk '{print substr($1, 1, length($1)-2)}') hz1_cnx6_cor${i}_scafftochr18ref.txt > hz1_cnx6_cor${i}.win.chr18.NONelevatedregion
grep -v $(awk -v ORS="\\\|" 1  elevatedregions.pos | awk '{print substr($1, 1, length($1)-2)}') hz1_cnx3P_cor${i}_scafftochr18ref.txt > hz1_cnx3P_cor${i}.win.chr18.NONelevatedregion
done

ENDTIME=$(date +%s)
echo $(date)
echo "It takes $(($ENDTIME - $STARTTIME)) seconds to complete this task"
echo "It takes $((($ENDTIME - $STARTTIME)/60)) mins to complete this task"
echo "It takes $((($ENDTIME - $STARTTIME)/3600)) hours to complete this task"
