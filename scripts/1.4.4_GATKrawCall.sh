#!/bin/bash -l
#SBATCH -J GATKrawCall
#SBATCH --get-user-env
#SBATCH --mail-user=gwee@biologie.uni-muenchen.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_production
#SBATCH --cpus-per-task=20
#SBATCH --time=14-00:00:00
#SBATCH -o /dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2/slurms/slurm-%j-%x.out

# sbatch 1.4.4_GATKrawCall.sh 00_134inds.txt 134inds GATKraw 2 chrDiploid poplist_134 DP3GQ0Miss10fullinfoQ100 3 0 0.1 100
#final selected filter: DP3GQ0Miss10fullinfoQ100, 134inds

module load vcftools/0.1.14-gcc8
module load bcftools/1.10.2-gcc8
conda activate biotools
#picard.jar=2.25.7

ref="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2/genome_HC_allpaths41687_v2.5_chrW.fasta"
dat="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2"
gatk="GATK_4.2.6.1"

echo $(date)
STARTTIME=$(date +%s)

cd ${dat}
cd ${dat}/05.1_recal/${3}

##### FILTER VCF
#fixed Info filter
QD=2
FS=60
SOR=3
MQ=40
MQRankSum=-12.5
ReadPosRankSum=-8
QUAL=100

#adjustable Geno filter
MINDPS=3
MINGQ=0
MINAD=3 #FMT/AD >= $MINAD

#F_Miss
MISS=0.1

java -Xmx16g -Djava.io.tmpdir=$dat/tmp -jar $dat/${gatk}.jar CombineGVCFs \
-R ${ref} \
$(cat $1) \
-O $2.gvcf.gz

#call including invariant (output required for dxy)
java -Xmx16g -Djava.io.tmpdir=$dat/tmp -jar $dat/${gatk}.jar GenotypeGVCFs \
-R ${ref} --max-alternate-alleles 4 \
-V $2.gvcf.gz -O $2_chrDiploid.vcf.gz -ploidy $4 --include-non-variant-sites 

#subset invariant sites
bcftools view --threads 20 -V indels --max-alleles 1 -Oz -o $2_$5.V1.gvcf.gz $2_chrDiploid.vcf.gz

#Geno level
bcftools filter --threads 16 -Oz -o $2_$5.V2.gvcf.gz --set-GTs . --include \
          "FMT/DP >= $MINDPS" \
          $2_$5.V1.gvcf.gz

bcftools view --threads 16 -Oz -o $2_$5.V3.gvcf.gz --include \
         "F_MISSING <= 0.25" \
         $2_$5.V2.gvcf.gz
bcftools index $2_$5.V3.gvcf.gz

# same filter on invariant and variant sites
##### FILTER VCF
QD=2
FS=60
SOR=3
MQ=40
MQRankSum=-12.5
ReadPosRankSum=-8
QUAL=100
MINDPS=3
MINGQ=0
MISS=0.1

#Info level
bcftools view --threads 20 -Oz -o $2_$5.V2_2.gvcf.gz --exclude \
       "QD < $QD || FS > $FS || SOR > $SOR || QUAL < $QUAL || MQ < $MQ || MQRankSum < $MQRankSum || ReadPosRankSum < $ReadPosRankSum" \
       $2_$5.V1.gvcf.gz

#Geno level
bcftools filter --threads 20 -Oz -o $2_$5.V3_2.gvcf.gz --set-GTs . --include \
         "FMT/DP >= $MINDPS" \
         $2_$5.V2_2.gvcf.gz

#remove missingness
bcftools view --threads 20 -Oz -o $2_$5.V4_2.gvcf.gz --include \
         "F_MISSING <= $MISS" \
         $2_$5.V3_2.gvcf.gz

# Remove repeat regions
vcftools --exclude-bed ${dat}/ref2.5repeats_header.bed --gzvcf $2_$5.V4_2.gvcf.gz --recode --out $2_$5_filtered_norepeats.gvcf
bcftools view -Oz -o $2_$5_filtered_norepeats.gvcf.gz $2_$5_filtered_norepeats.gvcf.recode.vcf
bcftools index $2_$5_filtered_norepeats.gvcf.gz
rm $2_$5_filtered_norepeats.gzvcf.recode.vcf

### generate variant sites only vcf ###
bcftools view --threads 20 -V indels --min-alleles 2 --max-alleles 2 -Oz -o $2_varonly_$5.V1.vcf.gz $2_$5_allsites.vcf.gz
bcftools view --threads 20 -V indels --min-alleles 3 -Oz -o $2_varonly_$5.multiallelic.vcf.gz $2_varonly_chrDiploid.vcf.gz
bcftools index $2_varonly_$5.multiallelic.vcf.gz
bcftools index -n $2_varonly_$5.multiallelic.vcf.gz

mkdir $2_$7
cd $2_$7

##### FILTER VCF
QD=2
FS=60
SOR=3
MQ=40
MQRankSum=-12.5
ReadPosRankSum=-8
QUAL=100
MINDPS=3
MINGQ=0
MINAD=3
MISS=0.1

#Info level
bcftools view --threads 20 -Oz -o $2_varonly_$5_$7.V2.vcf.gz --exclude \
       "QD < $QD || FS > $FS || SOR > $SOR || QUAL < $QUAL || MQ < $MQ || MQRankSum < $MQRankSum || ReadPosRankSum < $ReadPosRankSum" \
       ../$2_varonly_$5.V1.vcf.gz

#Geno level
bcftools filter --threads 20 -Oz -o $2_varonly_$5_$7.V3.vcf.gz --set-GTs . --include \
         "FMT/DP >= $MINDPS & FMT/GQ >= $MINGQ" \
         $2_varonly_$5_$7.V2.vcf.gz

#remove missingness
bcftools view --threads 20 -Oz -o $2_varonly_$5_$7_filtered.vcf.gz --include \
         "F_MISSING <= $MISS" \
         $2_varonly_$5_$7.V3.vcf.gz

echo "count filtered sites at INFO level of vcf"
SitesF=$(zcat $2_varonly_$5_$7.V2.vcf.gz | grep -v '#' | wc -l)
echo $SitesF
echo "count filtered sites at GENO level after missingness removed"
SitesG=$(zcat $2_varonly_$5_$7_filtered.vcf.gz | grep -v '#' | wc -l)
echo $SitesG

### Remove repeat regions
vcftools --exclude-bed ${dat}/ref2.5repeats_header.bed --gzvcf $2_varonly_$5_$7_filtered.vcf.gz --recode --out $2_varonly_$5_$7_filtered_norepeats
vcftools --vcf $2_varonly_$5_$7_filtered_norepeats.recode.vcf --singletons --out $2_varonly_$5_$7_filtered_norepeats
bcftools view -Oz -o $2_varonly_$5_$7_filtered_norepeats.vcf.gz $2_varonly_$5_$7_filtered_norepeats.recode.vcf
bcftools index $2_varonly_$5_$7_filtered_norepeats.vcf.gz
rm $2_varonly_$5_$7_filtered_norepeats.recode.vcf

#check mean quality scores
bcftools query -f '%CHROM %POS %REF %ALT %QUAL\n' *_norepeats.vcf.gz | awk '{ total += $5; count++ } END { print total/count }'

ENDTIME=$(date +%s)
echo $(date)
echo "It takes $(($ENDTIME - $STARTTIME)) seconds to complete this task"
echo "It takes $((($ENDTIME - $STARTTIME)/60)) mins to complete this task"
echo "It takes $((($ENDTIME - $STARTTIME)/3600)) hours to complete this task"
