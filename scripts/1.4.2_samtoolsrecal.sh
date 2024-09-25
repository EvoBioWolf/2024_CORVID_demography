#!/bin/bash -l
#SBATCH -J sam2recal
#SBATCH --get-user-env
#SBATCH --mail-user=gwee@biologie.uni-muenchen.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_production
#SBATCH --cpus-per-task=20
#SBATCH --time=8-00:00:00
#SBATCH -o /dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2/slurms/slurm-%j-%x.out

# sbatch 1.4.2_samtoolsrecal.sh 134inds samtools_DP3GQ0Miss10Q30 poplist_134_angsd 3 0 0.1 30
#final selected filter: samtools_DP3GQ0Miss10Q30 134inds

module load angsd/0.933-gcc8
#module load plink2

module load vcftools/0.1.14-gcc8
module load bcftools/1.10.2-gcc8
conda activate biotools
#picard.jar=2.25.7
#plink1.9 on conda

ref="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2/genome_HC_allpaths41687_v2.5_chrW.fasta"
dat="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2"

echo $(date)
STARTTIME=$(date +%s)

cd ${dat}/05.1_recal

minMapQ=30
minQ=30

parallel -j 30 '
minMapQ=30;
minQ=30;
file="bam_134inds_rescaled"
outdir="samtools/134inds/134inds"
ref="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2/genome_HC_allpaths41687_v2.5_chrW.fasta";
echo running scaffold_{};
bcftools mpileup -f ${ref} -b ${file}.filelist \
--redo-BAQ --annotate FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,FORMAT/SP,INFO/AD,INFO/ADF,INFO/ADR \
-q ${minMapQ} -Q ${minQ} -r scaffold_{} | bcftools call -mv -f GQ -Oz -o ./${outdir}_samtools_scaffold_{}.vcf.gz' ::: {0..1339}

cd ./samtools/$1
ls *scaffold_*.vcf.gz > vcf_list.txt
bcftools concat --file-list vcf_list.txt -Oz -o ${1}_samtools.vcf.gz --thread 30
bcftools view --threads 30 -V indels --min-alleles 2 --max-alleles 2 -Oz -o ${1}_samtools.V1.vcf.gz ${1}_samtools.vcf.gz

cd ${dat}/05.1_recal/samtools/$1
mkdir $2
cd $2

##### FILTER VCF
#fixed Info filter
MQ=40
MQRankSum=-12.5
ReadPosRankSum=-8
QUAL=30

#adjustable Geno filter
MINDPS=3
MINGQ=0

#F_Miss
MISS=0.1

#Info level
bcftools view --threads 20 -Oz -o ${1}_${2}.V2.vcf.gz --exclude \
        "QUAL < $QUAL || MQ < $MQ" \
       ../${1}_samtools.V1.vcf.gz

#Geno level
bcftools filter --threads 20 -Oz -o ${1}_${2}.V3.vcf.gz --set-GTs . --include \
          "FMT/DP >= $MINDPS & FMT/GQ >= $MINGQ" \
          ${1}_${2}.V2.vcf.gz

#remove missingness
bcftools view --threads 20 -Oz -o ${1}_${2}_filtered.vcf.gz --include \
         "F_MISSING <= $MISS" \
         $1_$2.V3.vcf.gz

bcftools index ${1}_${2}_filtered.vcf.gz

#sort before remove repeat region (to verify that the previous vcf is the same)
java -Xmx20g -Djava.io.tmpdir=$dat -jar $dat/picard.jar SortVcf TMP_DIR=$dat/05.1_recal/samtools \
      MAX_RECORDS_IN_RAM=10000 \
      I=${1}_${2}_filtered.vcf.gz \
      O=${1}_${2}_filtered_sorted.vcf
bgzip ${1}_${2}_filtered_sorted.vcf
vcftools --exclude-bed ${dat}/ref2.5repeats_header.bed --gzvcf ${1}_${2}_filtered_sorted.vcf.gz --recode --out $1_$2_filtered_norepeats_sorted
bcftools view -Oz -o $1_$2_filtered_norepeats_sorted.vcf.gz $1_$2_filtered_norepeats_sorted.recode.vcf
rm $1_$2_filtered_norepeats_sorted.recode.vcf
echo "unsorted_norepeats"
SitesG=$(zcat ${1}_${2}_filtered_norepeats.vcf.gz | grep -v '#' | wc -l)
echo $SitesG
echo "sorted_norepeats"
SitesS=$(zcat $1_$2_filtered_norepeats_sorted.vcf.gz | grep -v '#' | wc -l)
echo $SitesS

#Remove repeat regions
vcftools --exclude-bed ${dat}/ref2.5repeats_header.bed --gzvcf ${1}_${2}_filtered.vcf.gz --recode --out $1_$2_filtered_norepeats 
vcftools --vcf $1_$2_filtered_norepeats.recode.vcf --singletons --out $1_$2_filtered_norepeats
bcftools view -Oz -o $1_$2_filtered_norepeats.vcf.gz $1_$2_filtered_norepeats.recode.vcf
bcftools index $1_$2_filtered_norepeats.vcf.gz
rm $1_$2_filtered_norepeats.recode.vcf

#check mean quality scores
bcftools query -f '%CHROM %POS %REF %ALT %QUAL\n' *_norepeats.vcf.gz | awk '{ total += $5; count++ } END { print total/count }'

ENDTIME=$(date +%s)
echo $(date)
echo "It takes $(($ENDTIME - $STARTTIME)) seconds to complete this task"
echo "It takes $((($ENDTIME - $STARTTIME)/60)) mins to complete this task"
echo "It takes $((($ENDTIME - $STARTTIME)/3600)) hours to complete this task"
