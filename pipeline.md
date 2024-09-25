# Pipeline for the modern crow demography project
Path to folder containing all scripts: `./scripts`

# Raw reads processing
Raw reads were treated with Cutadapt 1.18 (M. Martin, 2011) to remove flanking adapter and reassessed with FastQC (Andrews, 2010) to ensure adapter contamination had been removed. 

`1.1.1_rawreads_pro.sh`
```
#!/bin/bash -l
#SBATCH -J cutadapt
#SBATCH --get-user-env
#SBATCH --mail-user=gwee@biologie.uni-muenchen.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_production
#SBATCH --cpus-per-task=8
#SBATCH --time=5-00:00:00
#SBATCH -o /dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/03_fresh/slurms/slurm-%j.out

# sbatch 1.1.1_rawreads_pro.sh (~5 days)
 
conda activate biotools 
#cutadapt=1.18
#fastqc=0.11.9 
#multiqc=0.9.1a0

dat="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2"

echo $(date)
STARTTIME=$(date +%s)

cd $dat

#trim adapter using cutadapt
cd rawreads_118
for fname in *pair_1*.fq.gz
do
base=${fname%.fq.gz*}
mv ${base}.fq.gz ${base}.R1.fq.gz
done

for fname in *pair_2*.fq.gz
do
base=${fname%.fq.gz*}
mv ${base}.fq.gz ${base}.R2.fq.gz
done

rename "pair_1_Illumina" "Illumina" *
rename "pair_2_Illumina" "Illumina" *

for fname in *.R1.fq.gz
do
base=${fname%.R1.fq.gz*}
cutadapt -j 8 -m 30 -q 20 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
-o ../01_trimmed/${base}_trimmed.R1.fq.gz -p ../01_trimmed/${base}_trimmed.R2.fq.gz ${base}.R1.fq.gz ${base}.R2.fq.gz
done

#do after all files have been trimmed
cd ../01_trimmed
mkdir 00_fastqc
for fname in *.fq.gz; do base=${fname%.fq.gz*}; fastqc ./${base}.fq.gz -t 8 --outdir ./00_fastqc; done
multiqc 00_fastqc
```

The trimmed reads were then mapped to the hooded crow reference genome (GenBank accession no. GCA_000738735.1) with BWA 0.7.17 (Li, 2013) and converted from SAM to BAM format using SAMtools 1.7 (Li et al., 2009). Read group information were added and dupli-cates were removed using picard 2.25.7 (http://broadinstitute.github.io/picard). 

`1.2.1_map.sh`
```
#!/bin/bash -l
#SBATCH -J map
#SBATCH --get-user-env
#SBATCH --mail-user=gwee@biologie.uni-muenchen.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=6
#SBATCH --time=2-00:00:00
#SBATCH -o /dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2/slurms/slurm-%j-%x.out

# for i in *.R2.fq.gz; do base=${i%_trimmed.R2.fq.gz*}; sbatch ../1.2.1_map.sh ${base}; done 

conda activate biotools 
#bwa=0.7.17
#samtools=1.7-1
#picard.jar=2.25.7

ref="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2/genome_HC_allpaths41687_v2.5_chrW.fasta"
dat="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2"

echo $(date)
STARTTIME=$(date +%s)
cd $dat

echo bwamem
cd 01_trimmed
bwa mem -t 6 -M ${ref} ${1}_trimmed.R1.fq.gz ${1}_trimmed.R2.fq.gz > ../02_mapped/${1}.sam 

echo samtools
cd ${dat}/02_mapped
samtools view -bhq 6 -@ 12 ${1}.sam > ${1}.bam
samtools sort -m 4G -@ 6 ${1}.bam -o ${1}.sorted.bam

mv ${1}.sorted.bam ./bam_sorted
mv ${1}.sam ./sam
mv ${1}.bam ./bam
```

`1.3.1_markdup.sh`
```
#!/bin/bash -l
#SBATCH -J markdup
#SBATCH --get-user-env
#SBATCH --mail-user=gwee@biologie.uni-muenchen.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=8
#SBATCH --time=2-00:00:00
#SBATCH -o /dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2/slurms/slurm-%j-%x.out

# cd 03_readgrp
# ls C.corone_*.bam | awk -F "_" '{print $1"_"$2}' | uniq > 00_C.corone.txt
# ls C.cornix_*.bam | awk -F "_" '{print $1"_"$2}' | uniq > 00_C.cornix.txt
# ls C.orientalis_*.bam | awk -F "_" '{print $1"_"$2}' | uniq > 00_C.orientalis.txt
# ls C.torquatus_*.bam | awk -F "_" '{print $1"_"$2}' | uniq > 00_C.torquatus.txt
# ls Hybrid_*.bam | awk -F "_" '{print $1"_"$2}' | uniq > 00_Hybrid.txt
# while read i; do sbatch ../1.3.1_markdup.sh $(echo $i); done < 00_C.corone.txt
# while read i; do sbatch ../1.3.1_markdup.sh $(echo $i); done < 00_C.cornix.txt
# while read i; do sbatch ../1.3.1_markdup.sh $(echo $i); done < 00_C.orientalis.txt
# while read i; do sbatch ../1.3.1_markdup.sh $(echo $i); done < 00_C.torquatus.txt
# while read i; do sbatch ../1.3.1_markdup.sh $(echo $i); done < 00_Hybrid.txt

# while read i; do sbatch ../1.3.1_markdup.sh $(echo $i); done < 00_C.torquatus_toepad.txt

conda activate biotools 
#bwa=0.7.17
#samtools=1.7-1
#picard.jar=2.25.7


ref="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2/genome_HC_allpaths41687_v2.5_chrW.fasta"
dat="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2"

echo $(date)
STARTTIME=$(date +%s)

cd $dat
cd $dat/03_readgrp

echo $1
ls $PWD/$1* | awk '{print "INPUT="$1}' > bam_list_$1

java -Djava.io.tmpdir=$dat -jar $dat/picard.jar MarkDuplicates \
MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 TMP_DIR=$dat \
$(cat bam_list_$1) \
REMOVE_DUPLICATES=TRUE \
OUTPUT=$dat/04_markdup/$1_markdup.bam \
METRICS_FILE=$dat/04_markdup/$1_markdup.txt

cd $dat/04_markdup
echo cleanreadsends
java -Djava.io.tmpdir=$dat -jar $dat/picard.jar CleanSam \
-I $1_markdup.bam -O $1_markdup_cleaned.bam -R ${ref} --CREATE_INDEX true

echo buildbamindex
java -Xmx4g -jar $dat/picard.jar BuildBamIndex INPUT=$1_markdup_cleaned.bam

#Summary Stats & Coverage
java -Xmx4g -jar $dat/picard.jar CollectAlignmentSummaryMetrics \
-R ${ref} \
-I $1_markdup_cleaned.bam \
-O ./00_stats/${1}.alignments.txt \
--METRIC_ACCUMULATION_LEVEL=SAMPLE \
--METRIC_ACCUMULATION_LEVEL=READ_GROUP

#calculate coverage
mosdepth --threads 8 --use-median --by 100000 --fast-mode --no-per-base ./00_stats/${1} ${1}_markdup_cleaned.bam
```
## Damage pattern of pectoralis 
DNA damage was assessed with MapDamage2.0 (Jónsson et al., 2013) and was detected in the toe pad samples with elevated levels of C to T and G to A deamination for the 5’ and 3’ ends of reads, thus the two toepad samples were further trimmed with Cutadapt to remove 5bp from the beginning of each forward and reverse reads. One of the toepad samples was ultimate-ly removed as deamination persisted when assessed with mapDamage 2.0. 

`1.3.2_mapdamage.sh`
```
#!/bin/bash -l
#SBATCH -J mapdam
#SBATCH --get-user-env
#SBATCH --mail-user=gwee@biologie.uni-muenchen.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=2
#SBATCH --time=2-00:00:00
#SBATCH -o /dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2/slurms/slurm-%j-%x.out

#cd 04*
#for i in C.torquatus_X*_cleaned.bam; do base=${i%.bam*}; sbatch ../1.3.2_mapdamage.sh ${base}; done

conda activate mapdamage

ref="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2/genome_HC_allpaths41687_v2.5_chrW.fasta"
dat="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2"

echo $(date)
STARTTIME=$(date +%s)

cd $dat/04_markdup
mapDamage -i $1.bam -r ${ref}
mapDamage -i $1.bam -r ${ref} -d results_$1 --plot-only

mv ./results_$1/$1.rescaled.bam ./ 
java -Xmx4g -jar $dat/picard.jar BuildBamIndex INPUT=$1.rescaled.bam
```

Sort by ploidy 

`1.3.4_ploidy.sh`
```
#!/bin/bash -l
#SBATCH -J ploidy
#SBATCH --get-user-env
#SBATCH --mail-user=gwee@biologie.uni-muenchen.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=8
#SBATCH --time=2-00:00:00
#SBATCH -o /dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2/slurms/slurm-%j-%x.out

# cd 04_markdup
# for fname in *_markdup_cleaned.bam; do base=${fname%.bam*}; sbatch ../1.3.4_ploidy.sh ${base}; done
# for fname in C.torquatus_X*_markdup_cleaned.bam; do base=${fname%.bam*}; sbatch ../1.3.4_ploidy.sh ${base}; done
# for fname in C.torquatus_X*.rescaled.bam; do base=${fname%.bam*}; sbatch ../1.3.4_ploidy.sh ${base}; done

conda activate biotools 
module load vcftools/0.1.14-gcc8
module load bcftools/1.10.2-gcc8
#picard.jar=2.25.7
#GATK=4.2.1

ref="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2/genome_HC_allpaths41687_v2.5_chrW.fasta"
dat="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2"
gatk="GATK_4.2.6.1"
kvcf="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2/knownsites_filtered.V2.vcf.gz"

echo $(date)
STARTTIME=$(date +%s)

cd ${dat}/04_markdup

#extract chrM, chrZ and chrW from bam
#chrZ of the 100 males should be ploidy2 but for now lumped with ploidy1 
samtools view -h $1.bam -@ 8 $(cat ${dat}/chrHaploid.scaffolds.list) -o $1.ploidy1.bam
samtools view -h $1.bam -@ 8 $(cat ${dat}/chrDiploid.scaffolds.list) -o $1.ploidy2.bam
samtools index -@ 8 $1.ploidy1.bam
samtools index -@ 8 $1.ploidy2.bam
```
# Variant calling
Variant calling was conducted using three callers: GATK 4.2.6.1 (McKenna et al., 2010); BCFtools 1.10.2 (Danecek et al., 2021); and ANGSD 0.933 (Korneliussen et al., 2014). We found 82.1% of the variants were identified by all three callers (Figure S1), and a total of 14,935,463 single nucleotide polymorphisms (SNPs) remained after low-quality variants, indels and repeated regions were removed.

## ANGSD
`1.4.1_angsdrecal.sh`
```
# sbatch 1.4.1_angsdrecal.sh 134inds_rescaled angsdrecalq30mindepth400genodepth3 poplist_134_angsd 121 400

#final selected filter: angsdrecalq30mindepth400genodepth3 134inds

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

cd $dat/05.1_recal

MinDepth=400
MaxDepth=6800
MinDepthInd=3
MinInd=121
minMapQ=30
minQ=30
genoMinDepth=3

angsd -out ./angsd/$1_$2_min$4 -bam bam_$1.filelist -GL 2 -doMaf 2 -doMajorMinor 1 -nThreads 30 -doCounts 1 -doDepth 1 -SNP_pval 1e-6 -ref ${ref} \
-doPlink 2 -doGeno -4 -doPost 1 -geno_minDepth ${genoMinDepth} \
-setMinDepth ${MinDepth} -setMaxDepth ${MaxDepth} -setMinDepthInd ${MinDepthInd} -minInd ${MinInd} \
-minMapQ ${minMapQ} -minQ ${minQ} 

cd angsd 
plink --tfile $1_$2_min$4 --allow-no-sex --allow-extra-chr --make-bed --out $1_$2_min$4
plink --bfile $1_$2_min$4 --allow-extra-chr --recode vcf --out $1_$2_min$4
plink --bfile $1_$2_min$4 --allow-extra-chr --missing
bgzip $1_$2_min$4.vcf
bcftools index $1_$2_min$4.vcf.gz

### Remove repeat regions
vcftools --exclude-bed ${dat}/ref2.5repeats_header.bed --gzvcf $1_$2_min$4.vcf.gz --recode --out $1_$2_min$4_norepeats
vcftools --vcf $1_$2_min$4_norepeats.recode.vcf --singletons --out $1_$2_min$4_norepeats
bcftools view -Oz -o $1_$2_min$4_norepeats.vcf.gz $1_$2_min$4_norepeats.recode.vcf
bcftools index $1_$2_min$4_norepeats.vcf.gz
rm $1_$2_min$4_norepeats.recode.vcf
```

## SAMtools
`1.4.2_samtoolsrecal.sh`
```
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
outdir="samtools/134inds"
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
```

## GATK

`1.4.3_GATKrawHap.sh`
```
# for fname in  C.torquatus_X*_markdup_cleaned.rescaled.bam; do base=${fname%_markdup_cleaned.rescaled.bam*}; sbatch ../1.4.3_GATKrawHap.sh ${base}; done

conda activate biotools 
module load vcftools/0.1.14-gcc8
module load bcftools/1.10.2-gcc8
#picard.jar=2.25.7
#GATK=4.2.1

ref="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2/genome_HC_allpaths41687_v2.5_chrW.fasta"
dat="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2"
gatk="GATK_4.2.6.1"
kvcf="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2/knownsites_filtered.V2.vcf.gz"

echo $(date)
STARTTIME=$(date +%s)

cd ${dat}/04_markdup

echo HaplotypeCaller
java -Xmx16g -Djava.io.tmpdir=$dat/tmp -jar $dat/${gatk}.jar HaplotypeCaller \
--native-pair-hmm-threads 4 \
-I $1_markdup_cleaned.rescaled.ploidy2.bam \
-R $ref \
-O $dat/05.1_recal/GATKraw/$1.gvcf.gz -ERC GVCF \
-ploidy 2 --output-mode EMIT_ALL_CONFIDENT_SITES -mbq 30

java -Xmx16g -Djava.io.tmpdir=$dat/tmp -jar $dat/${gatk}.jar HaplotypeCaller \
--native-pair-hmm-threads 4 \
-I $1_markdup_cleaned.rescaled.ploidy1.bam \
-R $ref \
-O $dat/05.1_recal/GATKraw/ploidy1/$1.gvcf.gz -ERC GVCF \
-ploidy 1 --output-mode EMIT_ALL_CONFIDENT_SITES -mbq 30

java -Xmx16g -Djava.io.tmpdir=$dat/tmp -jar $dat/${gatk}.jar HaplotypeCaller \
--native-pair-hmm-threads 4 \
-I $1_markdup_cleaned.ploidy2.bam \
-R $ref \
-O $dat/05.1_recal/GATKraw/$1.gvcf.gz -ERC GVCF \
-ploidy 2 --output-mode EMIT_ALL_CONFIDENT_SITES -mbq 30

java -Xmx16g -Djava.io.tmpdir=$dat/tmp -jar $dat/${gatk}.jar HaplotypeCaller \
--native-pair-hmm-threads 4 \
-I $1_markdup_cleaned.ploidy1.bam \
-R $ref \
-O $dat/05.1_recal/GATKraw/ploidy1/$1.gvcf.gz -ERC GVCF \
-ploidy 1 --output-mode EMIT_ALL_CONFIDENT_SITES -mbq 30
```
Emit both variant and invariant sites for all samples, then separate the invariant (gzvcf) and variant sites (vcf). Apply the same filter for both files. The invariant gzvcf file will be used to generate neutral invariant sites (see below).   

`1.4.4_GATKrawCall.sh`
```
#final selected filter: DP3GQ0Miss10fullinfoQ100, 134inds
# sbatch 1.4.4_GATKrawCall.sh 00_134inds.txt 134inds GATKraw 2 chrDiploid poplist_134 DP3GQ0Miss10fullinfoQ100 3 0 0.1 100

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
134inds_varonly_chrDiploid.vcf.gz

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

### Remove repeat regions
vcftools --exclude-bed ${dat}/ref2.5repeats_header.bed --gzvcf $2_varonly_$5_$7_filtered.vcf.gz --recode --out $2_varonly_$5_$7_filtered_norepeats
vcftools --vcf $2_varonly_$5_$7_filtered_norepeats.recode.vcf --singletons --out $2_varonly_$5_$7_filtered_norepeats
bcftools view -Oz -o $2_varonly_$5_$7_filtered_norepeats.vcf.gz $2_varonly_$5_$7_filtered_norepeats.recode.vcf
bcftools index $2_varonly_$5_$7_filtered_norepeats.vcf.gz
rm $2_varonly_$5_$7_filtered_norepeats.recode.vcf
```
## Combine overlapping SNPs

`1.4.5_overlap.sh`
```
# sbatch 1.4.5_overlap.sh 134inds overlapped poplist_134

module load vcftools/0.1.14-gcc8
module load bcftools/1.10.2-gcc8
conda activate biotools

dat="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2"
gatkset="134inds_DP3GQ0Miss10fullinfoQ100"
samtoolsset="samtools_DP3GQ0Miss10Q30"
angsdset="134inds_rescaled_angsdrecalq30mindepth400genodepth3_min121"

cd ${dat}/05.1_recal/GATKraw/${gatkset}
echo "GATK_134inds"
bcftools query -f '%CHROM %POS\n' *_norepeats.vcf.gz | awk '{print $1"_"$2}' > ${dat}/05.1_recal/overlap/${gatkset}_pos.txt
	
cd ${dat}/05.1_recal/samtools/134inds/${samtoolsset}
echo "samtools_134inds"
bcftools query -f '%CHROM %POS\n' *_norepeats.vcf.gz | awk '{print $1"_"$2}' > ${dat}/05.1_recal/overlap/${samtoolsset}_pos.txt

cd ${dat}/05.1_recal/angsd
echo "angsds_134inds"
bcftools query -f '%CHROM %POS\n' ${angsdset}_norepeats.vcf.gz | awk '{print $1"_"$2}' > ${dat}/05.1_recal/overlap/${angsdset}_pos.txt

cd ${dat}/05.1_recal/overlap
conda activate renv
Rscript ${dat}/overlap.R

awk 'FS="_" {print $1"_"$2, $3}' overlapping_variants_pos.txt | awk 'NR>1' > overlap_allcallers.pos
vcftools --positions overlap_allcallers.pos --gzvcf ${dat}/05.1_*/samtools/${1}/${samtoolsset}/*_norepeats.vcf.gz --recode --out $1_$2 
vcftools --vcf $1_$2.recode.vcf --singletons --out $1_$2
bcftools view -Oz -o $1_$2.vcf.gz $1_$2.recode.vcf
rm $1_$2.recode.vcf

#!!! Note that ${samtoolsset} and the overlappedVCF is not sorted according to index and this is fixed after filtering by SortVCF!!!

##### FILTER VCF
MQ=40
QUAL=30
MINDPS=3
MISS=0.1

#Info level
bcftools view --threads 20 -Oz -o ${1}_${2}.V1.vcf.gz --exclude \
        "QUAL < $QUAL || MQ < $MQ" \
        $1_$2.vcf.gz

#Geno level
bcftools filter --threads 20 -Oz -o ${1}_${2}.V2.vcf.gz --set-GTs . --include \
          "FMT/DP >= $MINDPS" \
          ${1}_${2}.V1.vcf.gz

#remove missingness
bcftools view --threads 20 -Oz -o ${1}_${2}_filtered_norepeats.vcf.gz --include \
          "F_MISSING <= $MISS" \
          $1_$2.V2.vcf.gz

#annotate and sort + poplist also sorted alphabetically
bcftools annotate --set-id '%CHROM\_%POS\_%REF\_%FIRST_ALT' $1_$2_filtered_norepeats.vcf.gz -Oz -o $1_$2_filtered_norepeats_ID.vcf.gz
java -Xmx20g -Djava.io.tmpdir=$dat -jar $dat/picard.jar SortVcf TMP_DIR=$dat \
      MAX_RECORDS_IN_RAM=10000 \
      I=$1_$2_filtered_norepeats_ID.vcf.gz \
      O=$1_$2_filtered_norepeats_ID_sorted.vcf
rm $1_$2.V*.vcf.gz $1_$2_filtered_norepeats_ID.vcf.gz
mv $1_$2_filtered_norepeats_ID_sorted.vcf $1_$2_filtered_norepeats.vcf
bgzip $1_$2_filtered_norepeats.vcf
bcftools index $1_$2_filtered_norepeats.vcf.gz
vcftools --gzvcf $1_$2_filtered_norepeats.vcf.gz --freq --out  $1_$2_filtered_norepeats
awk 'NR>1 {print $1, $2}' OFS="\t"  $1_$2_filtered_norepeats.frq >  $1_$2_filtered_norepeats.pos
```
## Paralog
The site frequency spectrum (SFS) of 134 individuals was examined and a slight elevation at frequency 0.5 was observed, suggesting presence of paralogs. Loci that significantly differed from Hardy Weinberg equilibrium (p-value < 0.05) were identified with VCFtools 0.1.14 (Danecek et al., 2011) and removed, resulting in a final set of 14,822,221 SNPs for down-stream analysis.   
`1.4.5_paralogs.sh`
```
#!/bin/bash -l
#SBATCH -J paralogs
#SBATCH --get-user-env
#SBATCH --mail-user=gwee@biologie.uni-muenchen.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=2
#SBATCH --time=2-00:00:00
#SBATCH -o /dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2/slurms/slurm-%j-%x.out

# sbatch 1.4.5_paralogs.sh 134inds overlapped poplist_134

module load vcftools/0.1.14-gcc8
module load bcftools/1.10.2-gcc8
conda activate biotools

dat="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2"
gatkset="134inds_DP3GQ0Miss10fullinfoQ100"
samtoolsset="samtools_DP3GQ0Miss10Q30"
angsdset="134inds_rescaled_angsdrecalq30mindepth400genodepth3_min121"

cd ${dat}/05.1_recal/overlap

#identify paralogs using cnx3 population
vcftools --gzvcf $1_$2_filtered_norepeats.vcf.gz  --keep ${dat}/06_results/neutral/cnx3.pop --recode --out $1_$2_filtered_norepeats_cnx3
vcftools --vcf $1_$2_filtered_norepeats_cnx3.recode.vcf --hardy --out $1_$2_filtered_norepeats_cnx3
awk 'NR>1 {print $1, $2, $6<0.05}' $1_$2_filtered_norepeats_cnx3.hwe | grep -w 1 | awk '{print $1, $2}' > $1_$2_filtered_norepeats_cnx3_hwedis.pos
grep -v "scaffold_1026\|scaffold_1056\|scaffold_107\|scaffold_1223\|scaffold_246\|scaffold_261\|scaffold_271\|scaffold_305\|scaffold_320\|scaffold_373\|scaffold_458\|scaffold_60\|scaffold_70\|scaffold_750\|scaffold_78\|scaffold_927\|scaffold_971\|scaffold_995" $1_$2_filtered_norepeats_cnx3_hwedis.pos > $1_$2_filtered_norepeats_cnx3_hwedis_nochr18.pos
vcftools --gzvcf $1_$2_filtered_norepeats.vcf.gz --exclude-positions $1_$2_filtered_norepeats_cnx3_hwedis_nochr18.pos --recode --out $1_$2_filtered_norepeats_hwe
bcftools view -Oz -o $1_$2_filtered_norepeats_hwe.vcf.gz $1_$2_filtered_norepeats_hwe.recode.vcf
```
# Combine invariant sites and the final set of variant sites  
run `1.4.5_gvcf_combineinvariant.sh`
```
#!/bin/bash -l
#SBATCH -J gvcf
#SBATCH --get-user-env
#SBATCH --mail-user=gwee@biologie.uni-muenchen.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_production
#SBATCH --cpus-per-task=16
#SBATCH --time=10-00:00:00
#SBATCH -o /dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2/slurms/slurm-%j-%x.out

# sbatch 1.4.5_gvcf_combineinvariant.sh 134inds overlapped poplist_134
# combine invariant sites from GATK with the final variant vcf

module load vcftools/0.1.14-gcc8
module load bcftools/1.10.2-gcc8
conda activate biotools

dat="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2"
gatkset="134inds_DP3GQ0Miss10fullinfoQ100"
samtoolsset="samtools_DP3GQ0Miss10Q30"
angsdset="134inds_rescaled_angsdrecalq30mindepth400genodepth3_min121"

### Merge the filtered vcf and gvcf back together
cd ${dat}/05.1_recal/overlap
bcftools concat --threads 16 --allow-overlaps $1_$2_filtered_norepeats_hwe.vcf.gz ${dat}/05.1_recal/GATKraw/134inds_chrDiploid.V3.gvcf.gz -Oz -o $1_$2_filtered_norepeats_hwe_incinvariants.vcf.gz
bcftools index $1_$2_filtered_norepeats_hwe_incinvariants.vcf.gz
```

# PCA

Variants extracted by SNPRelate on `pca.R`

```
library(devtools)
library(ggplot2)
library(SNPRelate)
library(openxlsx)
library(dplyr)
library(psych)
library(ggpubr)

args <- commandArgs(trailingOnly = TRUE)

setwd("/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2")
poplist <- read.delim(file="poplist_134.txt")

setwd(args[2])

# PCA with SNPRelate #
# mean FST with SNPRelate #
colorBlindBlack8  <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
safe_colorblind_palette <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888", "#000000", "#E69F00", "#56B4E9")

vcf.fn <- (file=paste(args[1],".vcf.gz",sep=""))
snpgdsVCF2GDS(vcf.fn,"vcf.gds", method="copy.num.of.ref")
snpgdsSummary("vcf.gds")
vcf.gdsfile <- snpgdsOpen("vcf.gds")

# Fst estimation
v <- snpgdsFst(vcf.gdsfile, population=poplist$pop, method="W&C84", autosome.only=FALSE)
v$Fst
v$MeanFst
summary(v$FstSNP)

v2 <- snpgdsFst(vcf.gdsfile, population=poplist$pop, method="W&H02", autosome.only=FALSE)
v2$Fst
v2$MeanFst
v2$Beta
summary(v2$FstSNP)

vcf.pca <- snpgdsPCA(vcf.gdsfile, autosome.only=FALSE, missing.rate=NaN) ##use autosome.only=false to include all loci
names(vcf.pca)
#variance proportion (%)
pc.percent <- vcf.pca$varprop*100
head(round(pc.percent, 2))
print(pc.percent)
tab <- data.frame(sample.id = vcf.pca$sample.id,
EV1 = vcf.pca$eigenvect[,1], # the first eigenvector
EV2 = vcf.pca$eigenvect[,2], # the second eigenvector
EV3 = vcf.pca$eigenvect[,3], 
EV4 = vcf.pca$eigenvect[,4], 
EV5 = vcf.pca$eigenvect[,5], 
EV6 = vcf.pca$eigenvect[,6], 
EV7 = vcf.pca$eigenvect[,7], 
EV8 = vcf.pca$eigenvect[,8], 
EV9 = vcf.pca$eigenvect[,9], 
EV10 = vcf.pca$eigenvect[,10], 
stringsAsFactors = FALSE)
#extract vectors out for external plotting
tab$POP <- poplist$pop

setwd("/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2")
setwd(args[3])
write.csv(tab, file=paste(args[1],"snprelate.txt",sep="_"))

pdf(file=paste(args[1],"pc1vs2_pca.pdf",sep="_"))
ggplot(tab, group=POP) +
  geom_point(aes(x=EV1, y=EV2, color=POP, shape=POP), size=4) +
  scale_shape_manual(name="Population", values=c(15, 16, 17, 18, 3, 4, 9, 10, 11, 12, 13, 5, 15, 16, 17)) +
  scale_color_manual(name="Population", values=safe_colorblind_palette) +
  labs(x=paste("PC1 ",round(pc.percent[1],2),"%",sep=""), y=paste("PC2 ",round(pc.percent[2],2),"%",sep="")) +
  guides(color=guide_legend("Population"),fill=guide_legend("Population")) +
  ggtitle(args[1]) +
  theme_bw(base_size=12) +
  theme(axis.line = element_line(colour = "black"),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.border = element_blank(),
  panel.background = element_blank(),
  axis.title.x = element_text(size=15),
  axis.title.y = element_text(size=15),
  axis.text.x=element_text (size=12),
  axis.text.y=element_text (size=12),
  legend.text=element_text(size=12),
  legend.title=element_text(size=12)) 
dev.off()
```

# Admixture

`1.5.1_admixture.sh`
```
# for i in {1..10}; do sbatch 1.5.1_admixture.sh 134inds_overlapped_filtered_norepeats_ldpruned $i poplist_134 overlap; done

conda activate biotools
ref="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2/genome_HC_allpaths41687_v2.5_chrW.fasta"
dat="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2"

echo $(date)
STARTTIME=$(date +%s)

cd $dat/05.1_recal/$4
admixture --cv ${1}.bed $2 -j4 | tee log${2}.out 

cd ${dat}/06_results/admixture
mkdir $1
cd ./$1
mv $dat/05.1_recal/$4/log${2}.out ./
mv $dat/05.1_recal/$4/${1}.${2}.Q ./
mv $dat/05.1_recal/$4/${1}.${2}.P ./

conda activate renv
Rscript ${dat}/admixtureplot.R $1 $2 ${dat}/06_results/admixture/$1 ../$3.txt
```

# Admixtools
Compute f-statistics and infer admixture proportion using admixtools2 R-interface: https://uqrmaie1.github.io/admixtools/

Using <i> C. brachyrhynchos </i>(am) as an outgroup. Same results when a more distant outgroup, <i> C. moneduloides </i> (mon) was used instead. Added `poplist_134_AM.ind` for pop/ind info. I used vcf2eigenstrat.py to convert VCF to eigen format: https://github.com/mathii/gdc/blob/master/vcf2eigenstrat.py

Run `admixtools_1000gen.R` with `1.5.3_admixtools.sh`
```
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

# conda activate r4.2
cd ${dat}/06_results/admixtools
Rscript ${dat}/06_results/admixtools/admixtools_1000gen.R $2 $6 $4 $5
```
`admixtools_1000gen.R`
```
library(usethis)
library(devtools) 
library(igraph)
library(admixtools)
library(tidyverse)
library(magrittr)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)
setwd("/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2/06_results/admixtools")

genotype_data = args[1]
dir = args[1]
# extract_f2(args[1], args[1])
mypops = c("am","cor1","cor2","cor3","cnx1","cnx2","cnx3P","cnx3S","cnx4","cnx5","cnx6","ori1","ori2","ori3","pec1")
f2_blocks = f2_from_precomp(dir, afprod = TRUE, pops=mypops)

df <- data.frame(run_no=c(1), bestscore=c(2))
 for (i in args[3]:args[4]) {
  opt_results = find_graphs(f2_blocks, outpop = "am", numadmix = as.numeric(args[2]),
                                stop_gen = 1000, stop_gen2=50) 
  winner = opt_results %>% slice_min(score, with_ties = FALSE)
  winner$score[[1]]
  df[i,1] <- i 
  df[i,2] <- winner$score[[1]]
  saveRDS(winner, file = paste(args[1],"_allpop_adm",args[2],"_run_",i, sep=""),
          ascii = FALSE, version = NULL,
          compress = TRUE, refhook = NULL)

  file_name = paste(args[1],"_allpop_adm",args[2],"_run_",i,".tiff", sep="")
  tiff(file_name, units="in", width=5, height=5, res=300)
  print(plot_graph(winner$edges[[1]]))
  dev.off()
  }
 print(df)
 write.table(df, file=paste(args[1],"_allpop_adm",args[2],"_bs",args[3],"to",args[4],"_bestscores.txt",sep=""), quote=FALSE, row.names = FALSE, sep="\t")
 ```
 
# Dsuite
Identify admixture using ABBA-BABA approach with Dsuite: https://github.com/millanek/Dsuite

Added `SETS_fbranch.txt`, `SETS_trios2.txt`, and `tree.nwk` for Dtrios, Dinvestigate and fbranch analysis. 

Run `1.5.4_dsuite.sh`
```
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

cd $dat/06_results/dsuite
cd $dat/06_results/dsuite
./Dsuite/Build/Dsuite Dtrios ${dat}/${1}.vcf.gz SETS_$2.txt -t tree.nwk
./Dsuite/Build/Dsuite Dinvestigate ${dat}f/${1}.vcf.gz SETS_$2.txt SETS_trios2.txt -n $2
mv *_$2_50_25.txt ./dinvestigate

./Dsuite/Build/Dsuite Fbranch tree.nwk SETS_$2_tree.txt > fbranch_matrix.txt
#plotting require python3.8
conda activate py3.8
python ./Dsuite/utils/dtools.py fbranch_mfor i in *.pop; do base=${i%.pop*}; sbatch /dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2/1.5.5_summarystats.sh 134inds_overlapped_filtered_norepeats_hwe_incinvariants ${base}; doneatrix.txt tree.nwk

# cd dinvestigate
mv cnx6_cnx*_cor*_50_25.txt ./cnx6_cnx_cor
mv cnx*_cnx6_pec1_50_25.txt ./cnx_cnx6_pec1
mv cor*_cor*_cnx*_50_25.txt ./cor_cor_cnx
mv ori*_ori*_cnx*_50_25.txt ./ori_ori_cnx
mv ori*_ori*_pec1_50_25.txt ./ori_ori_pec1
```

# Summary statistics

## Tajima's D and pi with vcftools

Run `1.5.5_summarystats.sh`

```
#!/bin/bash -l
#SBATCH -J sumstats
#SBATCH --get-user-env
#SBATCH --mail-user=gwee@biologie.uni-muenchen.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=2
#SBATCH --time=2-00:00:00
#SBATCH -o /dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2/slurms/slurm-%j-%x.out

# for i in *.pop; do base=${i%.pop*}; sbatch /dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2/1.5.5_summarystats.sh 134inds_overlapped_filtered_norepeats_hwe_incinvariants ${base}; done


conda activate biotools
module load bcftools/1.10.2-gcc8
module load vcftools/0.1.14-gcc8
#picard.jar=2.25.7

ref="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2/genome_HC_allpaths41687_v2.5_chrW.fasta"
dat="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2"

echo $(date)
STARTTIME=$(date +%s)

cd $dat/06_results
mkdir sumstats
cd sumstats
mkdir inc_invariants_sites
cd inc_invariants_sites

#compute persite pi etc including invariant sites
vcftools --gzvcf ${dat}/05.1_recal/overlap/$1.vcf.gz --keep $2.pop --recode --stdout | gzip -c > $2_incinvariants.vcf.gz
vcftools --gzvcf $2_incinvariants.vcf.gz --site-pi --out $2
vcftools --gzvcf $2_incinvariants.vcf.gz --TajimaD 1 --out $2
vcftools --gzvcf $2_incinvariants.vcf.gz --hardy --out $2
vcftools --gzvcf $2_incinvariants.vcf.gz --het --out $2

#window-based over 50kb windows 
vcftools --gzvcf $2.vcf.gz --window-pi 50000 --out $2_50kb
vcftools --gzvcf $2.vcf.gz --TajimaD 50000 --out $2_50kb

vcftools --gzvcf $2.vcf.gz --window-pi 10000 --out $2_10kb
vcftools --gzvcf $2.vcf.gz --TajimaD 10000 --out $2_10kb
```
To get the mean stats once all populations are done:
```
for i in *10kb.Tajima.D; do  echo $i; awk '{ sum += $4; n++ } END { if (n > 0) print sum / n; }' ${i}; done
for i in *10kb.windowed.pi; do  echo $i; awk '{ sum += $5; n++ } END { if (n > 0) print sum / n; }' ${i}; done
```

| Stats | SPA | GER | FRA | ITA | BUL | POL | SWE | RUS | COR | IRQ | ori1 | ori2 | ori3 | pec1 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 50kb TajD | -0.0785 | -1.02 | -0.911 | -0.691 | -0.78 | -0.722 | -0.973 | -0.58 | 0.596 | -0.112 | -0.652 | 0.543 | -0.618 | -0.218 |
| 10kb TajD | -0.0727 | -0.983 | -0.884 | -0.668 | -0.75 | -0.689 | -0.935 | -0.558 | 0.558 | -0.112 | -0.627 | 0.514 | -0.599 | -0.213 |
| 50kb pi/e-3 | 3.83 | 5.61 | 4.7 | 4.73 | 4.4 |  4.72 | 5.25 | 3.45 | 1.64 | 2.54 | 3.69 | 1.84 | 3.8 | 2.21 |
| 10kb pi/e-3 | 3.94 | 5.77 | 4.84 | 4.87 | 4.53 |  4.86 | 5.41 | 3.55 | 1.7 | 2.62 | 3.8 | 1.9 | 3.91 | 2.27 |

## Fst, Dxy, pi using Simon Martin's popgenWindows.py

(clone this): https://github.com/simonhmartin/genomics_general

<b>!!! Note the VCF used to generate the geno input should also include invariant sites for dxy and pi to be correctly computed. </b> 
Run `1.5.6_fst.sh`

```
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

module load vcftools
module load bcftools
conda activate biotools

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
```

```
#calculate mean fst without plotting
##nan removed, negative values converted to 0 and unmapped scaffolds to chr are removed
for i in *win.chr
do
echo $i
grep chr ${i} | awk '{print $10}' | grep -v "nan" | awk '$1<0{$1=0}1' | awk '{_+=$1}END{printf "Avg: %0.4f\n",_/NR}'
done > 00_summary_fst_all.txt
```

Run `1.5.6_tajD.sh`
#doesnt work with the invariant files

```
#!/bin/bash -l
#SBATCH -J tajD
#SBATCH --get-user-env
#SBATCH --mail-user=gwee@biologie.uni-muenchen.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=2
#SBATCH --time=2-00:00:00
#SBATCH -o /dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2/slurms/slurm-%j-%x.out

# for i in 10000 50000
# do
# sbatch 1.5.6_tajD.sh 134inds_overlapped_filtered_norepeats_hwe_incinvariants poplist_134 $i cor1
# # done

module load vcftools
module load bcftools
conda activate biotools

ref="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2/genome_HC_allpaths41687_v2.5_chrW.fasta"
dat="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2"
scaff="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/01_probes/scaffold2chr"

echo $(date)
STARTTIME=$(date +%s)

cd $dat/06_results
mkdir tajD
cd $dat/06_results/tajD

#Simon Martin's script
python ${dat}/06_results/fst/genomics_general/popgenWindows.py --popsFile $dat/06_results/fst/${2}.txt -g $dat/06_results/fst/${1}.geno.gz \
-f phased -w ${3} -p ${4} -T 4 --ploidy 2 -o TajD_${3}_${4}.csv.gz --analysis popFreq
```
---

# Neutral sites
Generate neutral sites for demographic analysis using `scripts_rwilliamson/vcfSummarizer.py`

Run `1.6.1_neutral.sh`
```
#!/bin/bash -l
#SBATCH -J neutral
#SBATCH --get-user-env
#SBATCH --mail-user=gwee@biologie.uni-muenchen.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=2
#SBATCH --time=2-00:00:00
#SBATCH -o /dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2/slurms/slurm-%j-%x.out
#SBATCH --mem-per-cpu=4763mb

# for i in *.pop; do base=${i%.pop*}; sbatch /dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2/1.6.1_neutral.sh 134inds_overlapped_filtered_norepeats_hwe ${base}; done

conda activate py2
module load vcftools/0.1.14-gcc8

ref="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2/genome_HC_allpaths41687_v2.5_chrW.fasta"
dat="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2"

echo $(date)
STARTTIME=$(date +%s)

cd $dat/06_results
mkdir neutral
cd neutral

#keep only varaint sites of each pop
vcftools --gzvcf ${dat}/05.1_recal/overlap/$1.vcf.gz --keep $2.pop --mac 1 --recode --stdout | gzip -c > $2.vcf.gz
zcat $2.vcf.gz | grep -v "INFO=<ID=\|FORMAT=<ID" | gzip -c > $2_edited.vcf.gz
rm $2.vcf.gz
mv $2_edited.vcf.gz $2.vcf.gz

#Generate a vcf summary using vcfSummarizer.py 
python ./scripts_rwilliamson/vcfSummarizer.py ./$2.vcf.gz \
${neu}/genome_HC_allpaths41687_v2.5_annotated.sites -q 0 -d 0 -D 1000 -L 0 -N 0 > $2.summary

#identify synonymous sites on each pop vcf
awk 'NR<=12 || $8==3' $2.summary > ./0fold/${2}_0fold.txt
awk 'NR<=12 || $8==4' $2.summary > ./4fold/${2}_4fold.txt
awk 'NR<=12 || $8==0' $2.summary > ./intergene/${2}_intergene.txt
awk 'NR<=12 || $8==1' $2.summary > ./intron/${2}_intron.txt
awk 'NR<=12 || $8==2' $2.summary > ./exon/${2}_exon.txt

cd $dat/06_results/neutral/intron
awk 'NR>12 {print $1, $2}' ${2}_intron.txt > ${2}.pos
cd $dat/06_results/neutral/intergene
awk 'NR>12 {print $1, $2}' ${2}_intergene.txt > ${2}.pos
cd $dat/06_results/neutral/4fold
awk 'NR>12 {print $1, $2}' ${2}_4fold.txt > ${2}.pos
```
Once neutral sites have been identified for all populations, concatenate all sites into a single positon file (1-based) and generate neutral vcf.

See ` 1.6.1_neutralvcf.sh`

```
#!/bin/bash -l
#SBATCH -J neutralvcf
#SBATCH --get-user-env
#SBATCH --mail-user=gwee@biologie.uni-muenchen.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=2
#SBATCH --time=2-00:00:00
#SBATCH -o /dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2/slurms/slurm-%j-%x.out
#SBATCH --mem-per-cpu=4763mb

# sbatch 1.6.1_neutralvcf.sh 134inds_overlapped_filtered_norepeats_hwe neuall_incIRQ
# sbatch 1.6.1_neutralvcf.sh 134inds_overlapped_filtered_norepeats_hwe_amcrow_moneduloides_monedula_biallele neuall_incIRQ

conda activate py2
module load vcftools/0.1.14-gcc8
#picard.jar=2.25.7

ref="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2/genome_HC_allpaths41687_v2.5_chrW.fasta"
dat="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2"

echo $(date)
STARTTIME=$(date +%s)

cd $dat/06_results/neutral

#once all done, do once

if [ $(ls neuall_incIRQ.pos | wc -l) -eq 1 ]
then
    echo "neutral position file already exist"
else
    echo "generate neutral without chr18 positions"
    cd $dat/06_results/neutral/intron
    cat *.pos | sort | uniq > all_pos.txt
    cd $dat/06_results/neutral/intergene
    cat *.pos | sort | uniq > all_pos.txt
    cd $dat/06_results/neutral/4fold
    cat *.pos | sort | uniq > all_pos.txt
    cd $dat/06_results/neutral
    cat intron/all_pos.txt intergene/all_pos.txt 4fold/all_pos.txt | sort | uniq > neuall_incIRQ_pre.pos
    grep -v "scaffold_1026\|scaffold_1056\|scaffold_107\|scaffold_1223\|scaffold_246\|scaffold_261\|scaffold_271\|scaffold_305\|scaffold_320\|scaffold_373\|scaffold_458\|scaffold_60\|scaffold_70\|scaffold_750\|scaffold_78\|scaffold_927\|scaffold_971\|scaffold_995" neuall_incIRQ_pre.pos > neuall_incIRQ.pos
    mkdir ${dat}/05.1_recal/overlap/neuall_incIRQ
    cp neuall_incIRQ.pos ${dat}/05.1_recal/overlap/neuall_incIRQ/neuall_incIRQ.pos
fi

cd ${dat}/05.1_recal/overlap/neuall_incIRQ
vcftools --positions $2.pos --gzvcf ${dat}/05.1_recal/overlap/${1}.vcf.gz --recode --out $1_$2
bcftools view -Oz -o $1_$2.vcf.gz $1_$2.recode.vcf
bcftools index $1_$2.vcf.gz
rm $1_$2.recode.vcf
bcftools index -n $1_$2.vcf.gz

vcftools --gzvcf $1_$2.vcf.gz --TsTv-summary --out $1_$2

```
I repeated similar script for invariant sites for the scaling step of demographic simulations. See `1.6.1_neutral_invariantsites.sh`  

# MSMC2

A mappability mask per scaffold of the genome
`2.1.1_SNPable.sh`
```
#sbatch 2.0_SNPable.sh ref2.5

conda activate biotools 
#bwa=0.7.17
#samtools=1.7-1
#picard.jar=2.25.7
#bcftools/1.10.2-gcc8

ref="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2/genome_HC_allpaths41687_v2.5_chrW.fasta"
dat="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2"
echo $(date)
STARTTIME=$(date +%s)

# download SNPable sourcefile, move to msmc folder
# make to compile
# run from the folder

cd $dat
cd $dat/msmc/SNPable

#below is to identify non-uniquely mapping regions using Heng Li's method 
./script/splitfa ${ref} 35 | split -l 20000000
mv x* short_reads

# cd short_reads
for i in x*
do
bwa aln -R 1000000 -O 3 -E 3 ${ref} ${i} > ../bwa_aln/${i}.sai
bwa samse ${ref} ../bwa_aln/${i}.sai ${i} > ../bwa_aln/${i}.sam
gzip ../bwa_aln/${i}.sam
done

cd $dat/msmc/SNPable/bwa_aln
for fname in *.sam.gz
do
base=${fname%.sam.gz*}
gzip -dc ${base}.sam.gz | ../script/gen_raw_mask.pl > ../encoded_mask/rawMask_35_${base}.fa #encoded
../script/gen_mask -l 35 -r 0.5 ../encoded_mask/rawMask_35_${base}.fa > ../encoded_mask/mask_35_50_${base}.fa #encoded
done

# cd ${dat}/msmc/SNPable/encoded_mask
cat mask_35_50_*.fa | sed 's/ 35 0.500//g' > ../${1}.mask_35_50_encoded.fa
sed 's/ 35 0.500//g' ../${1}.mask_35_50_encoded.fa > ../${1}.mask_35_50_encoded_edited.fa
rm  ../${1}.mask_35_50_encoded.fa

#modify python script to specify correct path to the mask file and output dir
cd $dat/msmc
python2 ./msmc-tools/makeMappabilityMask.py
```

`2.1.2_msmc_ind.sh`
```
# sbatch 2.1.2_msmc_ind.sh C.cornix_S05 cnx3 29.0x
# sbatch 2.1.2_msmc_ind.sh C.corone_E14 cor1 26.7x
# sbatch 2.1.2_msmc_ind.sh C.corone_D06 cor2 25.4x
# sbatch 2.1.2_msmc_ind.sh C.cornix_A12 cnx4 19.9x
# sbatch 2.1.2_msmc_ind.sh C.orientalis_A18 ori1 19.9x
# sbatch 2.1.2_msmc_ind.sh C.cornix_IRQ641341 cnx6 22.0x
# sbatch 2.1.2_msmc_ind.sh C.corone_FPa01 cor3 17.2x
# sbatch 2.1.2_msmc_ind.sh C.torquatus_X02 pec1 16.5
# sbatch 2.1.2_msmc_ind.sh C.cornix_B04 cnx2 14.7x
# sbatch 2.1.2_msmc_ind.sh C.cornix_P11 cnx3 14.1x
# sbatch 2.1.2_msmc_ind.sh C.cornix_IT12 cnx1 15.7x
# sbatch 2.1.2_msmc_ind.sh C.corone_FPd01 cor3 14.8x

# pip install utils

conda activate biotools 
module load bcftools
#bwa=0.7.17
#samtools=1.7-1
#picard.jar=2.25.7
#bcftools/1.10.2-gcc8

ref="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2/genome_HC_allpaths41687_v2.5_chrW.fasta"
dat="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2"
scaff="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/01_probes/scaffold2chr"

echo $(date)
STARTTIME=$(date +%s)

cd $dat
#each chr takes about 2 hrs (28chr+2=30*2=60hrs ~5days for each sample)

#only macrochr are used:
cd $scaff
cat chr1A.* chr1.* chr2.* chr3.* chr4.* chr5.* > macrochr.scaffolds
sed -i "s/scaffold_//g" macrochr.scaffolds
mv macrochr.scaffolds $dat/msmc

cd $dat/msmc
mkdir ${1}
for i in $(cat macrochr.scaffolds)
for i in {0..100}
do 
samtools mpileup -q 20 -Q 20 -C 50 -u -r scaffold_${i} -f ${ref} ${dat}/04_markdup/${1}_markdup_cleaned.ploidy2.bam | \
bcftools call -c -V indels | \
./msmc-tools/bamCaller.py $(echo $(samtools depth -r scaffold_$i ${dat}/04_markdup/${1}_markdup_cleaned.ploidy2.bam | awk '{sum += $3} END {print sum / NR}')) ./${1}/$1_scaffold_${i}_mask.bed.gz | \
gzip -c > ./${1}/$1_scaffold_$i.vcf.gz
done

for i in {101..114}
do
samtools mpileup -q 20 -Q 20 -C 50 -u -r scaffold_${i} -f ${ref} ${dat}/04_markdup/${1}_markdup_cleaned.ploidy2.bam | \
bcftools call -c -V indels | \
./msmc-tools/bamCaller.py $(echo $(samtools depth -r scaffold_$i ${dat}/04_markdup/${1}_markdup_cleaned.ploidy2.bam | awk '{sum += $3} END {print sum / NR}')) ./${1}/$1_scaffold_${i}_mask.bed.gz | \
gzip -c > ./${1}/$1_scaffold_$i.vcf.gz
done

####create input files for msmc####
for i in $(cat macrochr.scaffolds)
for i in {0..100}
do
./msmc-tools/generate_multihetsep.py --mask=./${1}/${1}_scaffold_${i}_mask.bed.gz \
--mask=${dat}/msmc/ref2.5_masked/ref2.5_scaffold_${i}.mask.bed.gz \
./${1}/${1}_scaffold_${i}.vcf.gz > ./${1}/${1}_scaffold_${i}_input.txt
done

for i in {101..114}
do
./msmc-tools/generate_multihetsep.py --mask=./${1}/${1}_scaffold_${i}_mask.bed.gz \
--mask=${dat}/msmc/ref2.5_masked/ref2.5_scaffold_${i}.mask.bed.gz \
./${1}/${1}_scaffold_${i}.vcf.gz > ./${1}/${1}_scaffold_${i}_input.txt
done

#11 scaffs are chrZ (thus filesize=0) and 4 scaffs are chr18
cd ./${1}
mkdir chr18
for i in $(cat ${scaff}/chr18.scaffolds)
do
mv *$i* ./chr18
done

####generate msmc results####
cd ${dat}/msmc
./msmc-tools/getStats $(ls -l ./${1}/${1}_scaffold_*_input.txt | awk '{if ($5 != 0) print $9}') > ./${1}/${1}_getstats.txt

for i in {1..2}
do
./msmc${i} -o ./results/${1}_msmc${i}_output_32segs $(ls -l ./${1}/${1}_scaffold_*_input.txt | awk '{if ($5 != 0) print $9}')
./msmc${i} -p 1*2+15*1+1*2 -o ./results/${1}_msmc${i}_output_19segs $(ls -l ./${1}/${1}_scaffold_*_input.txt | awk '{if ($5 != 0) print $9}')
done

#default -p:1*2+25*1+1*2+1*3 (32segs)
#after multiple runs, 22-20 segments with 2 groups (1*20+1*2) worked the best for MSMC1. but for MSMC2, less segmentation is required
```

# Stairway plot

1dSFS per population generated with easySFS 
`2.2.1_1dsfs.sh`
```
# for i in *.poplist
# do
# base=${i%.poplist*}
# sbatch /dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2/2.2.1_1dsfs.sh 05.1_recal/overlap/neuall_incIRQ/134inds_overlapped_filtered_norepeats_hwe_neuall_incIRQ_nomissing_ldpruned ${base} $(wc -l $i | sed "s/$i//g")
# done

module load vcftools/0.1.14-gcc8
module load bcftools/1.10.2-gcc8
conda activate easySFS

ref="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2/genome_HC_allpaths41687_v2.5_chrW.fasta"
dat="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2"

echo $(date)
STARTTIME=$(date +%s)

cd ${dat}/stairway

nseq=$(echo $(( 2*$3 )))

#folded
# python3 ${dat}/fastsimcoal2/easySFS/easySFS.py -i ${dat}/${1}.vcf.gz -p ${2}.poplist -a -f -o $2_folded --prefix $2_folded --proj ${nseq}
```

`2.2.2_stairway_scaled.sh`
```
# for i in *.poplist
# do
# base=${i%.poplist*}
# sbatch /dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2/2.2.2_stairway.sh ${base}
# done

conda activate biotools 

ref="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2/genome_HC_allpaths41687_v2.5_chrW.fasta"
dat="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2"

echo $(date)
STARTTIME=$(date +%s)

cd $dat/stairway

#edit blueprint
cd stairway_scaled

java -cp stairway_plot_es Stairbuilder ${1}_folded.blueprint
bash ${1}_folded.blueprint.sh
bash ${1}_folded.blueprint.plot.sh
```
## How is scaling carried out?

I edited the "L" paramter to a total of 329,623,861 bp instead of 1.2 billion as scaling should be applied on ldpruned and neutral variant set - see the blueprint files in the stairway folder

No. of invariant sites: 625,134,484 (no chr18) \
No. of ldpruned neutral SNPs: 5,774,276 \
No. of neutral SNPs: 11,146,221 (no chr18)

Therefore scaling applied:
5,774,276/11,146,221*(625,134,484+11,146,221) = 329,623,861 

`2.2.2_stairway_rescaled.sh`
```
# for i in *.poplist
# do
# base=${i%.poplist*}
# sbatch /dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2/2.2.2_stairway_rescaled.sh ${base}
# done

conda activate biotools 

ref="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2/genome_HC_allpaths41687_v2.5_chrW.fasta"
dat="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2"

echo $(date)
STARTTIME=$(date +%s)


cd $dat/stairway

#edit blueprint
cd stairway_scaled

java -cp stairway_plot_es Stairbuilder ${1}_folded.blueprint
bash ${1}_folded.blueprint.sh
bash ${1}_folded.blueprint.plot.sh
```

# SMC++

`2.3.1_smcpp_neuperpopvcf.sh`
```
for i in c*.poplist; do base=${i%.poplist*}; sbatch /dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2/2.3.1_*.sh 134inds_overlapped_filtered_norepeats_hwe_neuall_incIRQ_ldpruned ${base}; done

conda activate singularity
#there's a cap of 200 container images per 6 hours across all my submitted jobs to singular!
cd ${dat}/smcpp/vcf
for pop in $(cat run.poplist)
do
cd ./${pop}
for i in {1..5}
do
for fname in chr${i}_${pop}_scaffolds.txt
do
chr=${fname%_${pop}_scaffolds.txt*}
for scaff in $(sed "s/_${pop}.vcf.gz//g" ${fname})
do
singularity run --bind $PWD:$HOME docker://terhorst/smcpp:latest vcf2smc -v ${chr}_${pop}.vcf.gz ./smc_output/${pop}_${chr}_${scaff}.smc.gz ${scaff} ${pop}:$(cat ${dat}/smcpp/vcf/${pop}.poplist | paste -s -d, -)
done
done
done
sleep 21600
done

cd ${dat}/smcpp/vcf
for pop in $(cat run.poplist)
do
cd ./${pop}/smc_output
singularity run --bind $PWD/:$HOME docker://terhorst/smcpp:latest estimate 3.18e-09 $(ls -l ${pop}_chr*.smc.gz | awk '{if ($5 > 1000) print $9}')
cd ${dat}/smcpp/vcf
done

##combine plots, run only after all pops are done
#cp *.json ../results/${3}/${1}.json
#cd ../results/${3}
#singularity run --bind $PWD/:$HOME docker://terhorst/smcpp:latest plot smcpp_all_${3}.pdf *.json -g 5.79 -c
```

# Fastsimcoal2
Using fastsimcoal2.7 to infer the demography of the C. corone complex and compare two main models - whether central European crows (CEU) diverged from Spain (SPA) or the hooded crows (ESE). 

---

## Generate SFS

<b> Generate folded SFS </b>: 

Run `1.6.3_easySFS_folded.sh` 

Generates unfolded and folded SFS from the ancestral state ascertained VCF with `easySFS.py`. A total of 50 individuals (15 Spain, 15 CEU, 15 ESE and 5 IRQ) were randomly selected to produce js-SFS.
```
#!/bin/bash -l
#SBATCH -J easySFS
#SBATCH --get-user-env
#SBATCH --mail-user=gwee@biologie.uni-muenchen.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=2-00:00:00
#SBATCH -o /dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2/slurms/slurm-%j-%x.out

# sbatch 1.6.3_easySFS.sh 05.1_recal/overlap/neuall_incIRQ/134inds_overlapped_filtered_norepeats_hwe_neuall_incIRQ cor1_cor2to3_cnx1to3_cnx6_50ind_folded_all_nochr18 30,30,30,10 4PopModel1


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

#unfolded SFS
python3 ${dat}/fastsimcoal2/easySFS/easySFS.py -i ${dat}/134inds_overlapped_filtered_norepeats_hwe_amcrow_moneduloides_monedula_biallele_neuall_incIRQ_outgroup_noamcrow2_AAknown_final.vcf.gz -p ${2}.pop -a -f --unfolded --proj ${3} -o $2 --prefix ${4}

#folded SFS
python3 ${dat}/fastsimcoal2/easySFS/easySFS.py -i ${dat}/134inds_overlapped_filtered_norepeats_hwe_neuall_incIRQ.vcf.gz -p ${2}.pop -a -f --proj ${3} -o $2 --prefix ${4}
```

<b>IMPORTANT</b>:
make sure invariant sites are corrected in the SFS matrix output. I modified column 0,0 of each observed 2d-SFS to 625,134,484 invariant sites. This is crucial as Fastsimcoal does scaling internally. 

<b> No. of invariant sites</b>: 625,134,484 (no chr18)\
<b>No of of SNPs for folded SFS</b>: 11,146,221 (no chr18)\

---

## Fastsimcoal simulation
Path to input and output files: `fastsimcoal2/4Pop/cor1_cor2to3_cnx1to3_cnx6_50ind_folded_all/fastsimcoal2`

See the input files required for Fastsimcoal to run: `4PopModel7x.est` `4PopModel7x.tpl` `4PopModel8x.est` `4PopModel8x.tpl`

<b> Estimate parameters </b>: \
Finally we are ready to estimate parameters with the observed SFS and input files.

Run `1.6.3_fastsimcoal_folded.sh` to generate 100 replicates for each model with folded SFS:
 ```
#!/bin/bash -l
#SBATCH -J fsimcoal
#SBATCH --get-user-env
#SBATCH --mail-user=gwee@biologie.uni-muenchen.de
#SBATCH --clusters=serial
#SBATCH --partition=serial_std
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --time=4-00:00:00
#SBATCH -o /dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2/slurms/slurm-%j-%x.out

# for bs in {1..100};do for i in {7..8}; do sbatch 1.6.3_fastsimcoal_folded.sh cor1_cor2to3_cnx1to3_cnx6_50ind_folded_all_nochr18 30,30,30,10 4PopModel${i}x ${bs}; done; do$

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

```
<b> Best estimate for each model </b>: \
Out of the 100 replicates, we want the best estimates for each model. We also want to compare the likehoods and AIC of each model to assess which is a better fit of the observed data.  

Run `1.6.4_plot.sh`:
```
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


# sbatch 1.6.4_plot.sh cor1_cor2to3_cnx1to3_cnx6_50ind_unfolded_all 4PopModel4x3 100 cor1 cor2to3 cnx1to3 cnx6
# sbatch 1.6.4_plot.sh cor1_cor2to3_cnx1to3_cnx6_50ind_unfolded_all 4PopModel3x3 100 cor1 cor2to3 cnx1to3 cnx6

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

#smallest difference between observed and expected likelihoods
cat run{1..100}/${2}/${2}.bestlhoods | grep -v "MaxObsLhood\|MaxEstLhood" | awk '{print NR,$(NF-1),$NF}' | awk '{print NR,$3-$2}' | sed -e "s/-9/NA/g" | grep -v "NA" | sort$
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
#these rscripts are downloaded from the fastsimcoal package
Rscript ${dat}/fastsimcoal2/calculateAIC.sh $2
Rscript ${dat}/fastsimcoal2/SFStools.R -t print2D -i $2 -z
Rscript ${dat}/fastsimcoal2/SFStools.R -t 2Dto1D -i $2 -z
Rscript ${dat}/fastsimcoal2/SFStools.R -t print1D -i $2 -z
Rscript ${dat}/fastsimcoal2/plotModel.R -p $2 -l ${5},${6},${7},${8}

#see which model with lowest AIC
zless ./*/*.AIC > AIC.txt
```

### Results 

I ran a few variation of models to test if CEU is from Spain or hooded crows. Here is a summary of the AIC of each model for the folded and unfolded SFS

| Model comparison	| AIC for CEU from Spain |	AIC for CEU from Hooded | Remarks |	
| --------| -------- | -------- | -------- |
| <b> Folded (no chr18)</b> |
| 1 vs 2 | 430,198,123 | 430,350,551 |  no GR |
| 3 vs 4 | 430,218,015 | 430,354,379 |  same as 3x4 and 4x3 |
| 7x vs 8x | 430,169,400 | 430,364,442 |  2 GR of diff rate |

<b> Folded (without chr18) </b>

| 99/100 rep | Model1 | Model 2 |
:---:|:---:|:---:
|  | <img src="./fastsimcoal/folded/nochr18/4PopModel1.jpeg" width="600" height=auto /> | <img src="./fastsimcoal/folded/nochr18/4PopModel2.jpeg" width="600" height=auto /> |
| 93/98 rep | Model3 | Model 4 |
|  | <img src="./fastsimcoal/folded/nochr18/4PopModel3.jpeg" width="600" height=auto /> | <img src="./fastsimcoal/folded/nochr18/4PopModel4.jpeg" width="600" height=auto /> |
| 100/100 rep | Model7x | Model 8x |
|  | <img src="./fastsimcoal/folded/nochr18/4PopModel7x.jpeg" width="600" height=auto /> | <img src="./fastsimcoal/folded/nochr18/4PopModel8x.jpeg" width="600" height=auto /> |


### Summary of best models from jaatha and fsimcoal runs - without chr18

13 sep 24: updated with corrected Jaatha estimates
log-likelihoods: -4508742 vs. -5774682 
AIC: 20,763,566.23 vs 26,591,112
| Para | folded Jaatha_Spain | folded fsimcoal_Spain | folded Jaatha_hooded | folded fsimcoal_hooded |
| --- | --- | --- | --- | --- |
| NANC1 | 9.41e+05 | 1.134345e+06 | 2.56e+05  | 7.864780e+05 |
| NANC2	| 9.41e+05 | 1.080252e+06 | 8.08e+04 | 1.074962e+06 |
| NANC3 | 8.78e+04 | 7.633200e+04 | 7.28e+04 | 6.994000e+04 |
| NPOP1	| 1.65e+04 | 1.001700e+04 | 2.10e+04 | 1.757000e+04 |
| NPOP2	| 7.50e+04 | 2.743970e+05 | 2.35e+05 | 1.101240e+06 |
| NPOP3	| 1.11e+05 | 4.668310e+05 | 1.74e+05 | 5.062650e+05 |
| NPOP4	| 3.77e+04 | 3.556700e+04 | 6.64e+04 | 5.333400e+04 |
| TDIV3	| 7.39e+04 | 8.699200e+04 | 1.24e+05 | 9.652700e+04 |
| TDIV2	| 1.84e+04 | 2.265400e+04 | 1.03e+05 | 6.224900e+04 |
| TDIV1	| 5.39e+03 | 5.962000e+03 | 1.03e+05 | 2.064400e+04 |
| TMRMG	| 5.29e+03 | 3.455000e+03 | 9.57e+03 | 8.861000e+03 |
| TEGR1 | 4.84e+03 | 5.360000e+02 | 6.56e+04 | 7.626000e+03 |
| TEGR2 | 3.08e+03 | 5.962000e+03 | 1.27e+04 | 1.365700e+04 |
| GR1 | -3.74e-05 | -2.349845948608230e-04 | -8.84e-07 | -2.567120337669507e-04 |
| GR2 | -8.74e-05 | -2.536942175178945e-04 | -1.88e-05 | -1.961217770466957e-04 |
| MIG01R |  9.52e-05 | 1.686905001543047e-04 | 1.09e-04 | 1.474836491352024e-04 |
| MIG10R | 9.70e-05 | 3.161484787677673e-04 | 1.84e-05 | 1.742554869641556e-05 |
| MIG12R | 2.42e-04 | 6.942012945324323e-04 | 7.21e-05 | 4.963224397438161e-05 |
| MIG21R | 9.32e-05 | 1.294168316901755e-04 | 5.02e-05 | 5.461327335068656e-05 |
| MIG23R | 8.82e-06  | 5.676149703025700e-06 | 8.82e-06 | 5.793928931137984e-06 |
| MIG32R| 4.41e-05 | 4.524627243921338e-05 | 2.39e-05 | 2.921071603081070e-05 |
| p19 | 1.09e-09 |  | 1.39e-09 |  |

---

## Check model fitness
This step is actually conducted in the previous step using `SFStools.R -t print2D`. However, the SFS generated by Fastsimcoal sums up to 1 and the observed SFS is scaled accordingly. I will simulate the best estimates by fastsimcoal and Jaatha to SFS to assess the fitness of both to the observed data.

First, we need to prepare a par file for each of the best model. Copy the maxL.par output of the best model and edit the bottom to specify simulation of DNA sequences. I will simulate 636281 of 1000bp unlinked loci to scale to my folded input data (11,146,221 SNPs+625,134,484 invariant sites)

See input files in folder `fastsimcoal/4Pop/cor1_cor2to3_cnx1to3_cnx6_50ind_folded_all_nochr18/fastsimcoal2/bestruns/modelfit`
```
//Number of independent loci [chromosome] 
6362807 0 
//Per chromosome: Number of contiguous linkage Block: a block is a set of contiguous loci
1
//per Block:data type, number of loci, per generation recombination and mutation rates and optional parameters
DNA 100 0 3.18e-9 0.5
```

Run `1.6.4_modelfit_folded.sh` for folded SFS
```
#!/bin/bash -l
#SBATCH -J modelfit
#SBATCH --get-user-env
#SBATCH --mail-user=gwee@biologie.uni-muenchen.de
#SBATCH --clusters=cm2_tiny
#SBATCH --partition=cm2_tiny
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=3-00:00:00
#SBATCH -o /dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2/slurms/slurm-%j-%x.out

# for i in {7..8}; do sbatch 1.6.4_modelfit_folded.sh 4PopModel${i}_fastsimcoal 4PopModel${i} fastsimcoal cor1_cor2to3_cnx1to3_cnx6_50ind_folded_all_nochr18; done
# for i in {3..4}; do sbatch 1.6.4_modelfit_folded.sh 4PopModel${i}x_jaatha 4PopModel${i}x jaatha; done

conda activate easySFS

ref="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2/genome_HC_allpaths41687_v2.5_chrW.fasta"
dat="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2"

echo $(date)
STARTTIME=$(date +%s)

cd ${dat}/fastsimcoal2/4Pop/${4}/fastsimcoal2
cd bestruns/modelfit

#infinite sites
${dat}/fastsimcoal2/fsc27 -i ${1}_infinite.par -j -m -s0 -x -I -q -n100 -c 2

#finite sites
${dat}/fastsimcoal2/fsc27 -i ${1}_finite.par -j -m -s0 -x -q -n100 -c 2

##legend##
#-j: output one simulated or bootstrapped SFS per file in a separate directory for easier analysis
#-m: Computes the site frequency spectrum (SFS) for th.pare minor alleles for each population sample and for all pairs of samples (joint SFS)
#-s0: output all SNPs in the DNA sequence
#-x: Does not generate Arlequin output file
#-I:Generates DNA mutations according to an infinite site (IS) mutation model
#-q: quiet mode
#-n100: 100 simulations

#remove the "no. of removed sites" text on finite SFS
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
Rscript ${dat}/fastsimcoal2/SFStools.R -t print1D -i $2 -z
done
cd ../
done
```

---
## Parametric bootstrap
Now we want to assess the confidence interval of each best model by simulating on the 100 bootstrapped SFS generated in the previous step. I ran this with folded 3x3 model only.

<b> !!! Finite SFS </b>


Run `1.6.5_simparabs_folded.sh` for folded SFS: 
```
#!/bin/bash -l
#SBATCH -J parabsfolded
#SBATCH --get-user-env
#SBATCH --mail-user=gwee@biologie.uni-muenchen.de
#SBATCH --clusters=serial
#SBATCH --partition=serial_std
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=4-00:00:00
#SBATCH -o /dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2/slurms/slurm-%j-%x.out

for i in {1..100}; do sbatch 1.6.5_simparabs_folded.sh 4PopModel7_fastsimcoal $i finite; done 

conda activate easySFS

ref="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2/genome_HC_allpaths41687_v2.5_chrW.fasta"
dat="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2"

echo $(date)
STARTTIME=$(date +%s)

cd ${dat}/fastsimcoal2/4Pop/cor1_cor2to3_cnx1to3_cnx6_50ind_folded_all_nochr18/fastsimcoal2
cd bestruns/modelfit
cd ${1}_${3}

cp ../${1}.tpl ./${1}_${3}_${2}/${1}_${3}.tpl
cp ../${1}.est ./${1}_${3}_${2}/${1}_${3}.est
cp ../${1}.pv ./${1}_${3}_${2}/${1}_${3}.pv

cd ${1}_${3}_${2}
echo ${1}_${3}_${2}
${dat}/fastsimcoal2/fsc27 -t ${1}_${3}.tpl -e ${1}_${3}.est –initvalues ${1}_${3}.pv -d -M -L 100 -n 1000000 -q -c 8

```
---

### Bias correction
Run on R `boxplot_bscomparison.R`
```
library(ggplot2)

com <- c("Model7x")
header <- c("RUN","NPOP1","NPOP2","NPOP3","NPOP4","NANC1","NANC2","NANC3","TMRMG",
           "MIG01R","MIG10R" ,"MIG12R" ,"MIG21R","MIG23R","MIG32R","GROW1", "GROW2", "TDIV3","TDIV2","TDIV1", "TEGR1", "TEGR2", "MaxEstLhood","MaxObsLhood")

setwd("/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2/fastsimcoal2/4Pop/cor1_cor2to3_cnx1to3_cnx6_50ind_folded_all_nochr18/fastsimcoal2/bestruns")
standardmodelSPA <- read.delim(paste("./4Pop",com[1],"/4Pop",com[1],".bestlhoods",sep=""), header=TRUE, sep="\t")
RUN <- 0
standardmodelSPA <- cbind(RUN,standardmodelSPA)
setwd("/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2/fastsimcoal2/4Pop/cor1_cor2to3_cnx1to3_cnx6_50ind_folded_all_nochr18/fastsimcoal2/bestruns/modelfit/4PopModel7x_fastsimcoal_finite")
modelSPAbs <- read.delim(paste("4Pop",com[1],"_fastsimcoal_finite.bestlhoods.all",sep=""), header=FALSE, sep="\t")
colnames(modelSPAbs) <- header
modelSPA <- rbind(standardmodelSPA,modelSPAbs)
modelSPA <- modelSPA[,c(-23, -24)]

# bias correction

bc <-  function(t0, tt){2*t0-colMeans(tt)}
ci <- function(tt){quantile(tt, c(0.025, (1 - 0.025)))}
t0 = modelSPA[1,][2:(length(modelSPA))]
tt = modelSPAbs[2:(length(modelSPAbs)-2)]

for (i in 1:length(t0)) {
  print(bc(t0[i],tt[i]))
  print(ci(tt[,i]))  
}

# boxplot

pdf(file="model7xparabs_boxplot.pdf")
modelSPA$RUN <- c(TRUE, rep(FALSE, nrow(modelSPA) - 1L))  # tag first row as interesting
df.2 <- melt(modelSPA)  # convert df to long format
ggplot(subset(df.2, !RUN), aes(x=variable, y=value)) + 
  geom_boxplot() + scale_y_log10() +
  geom_point(data=subset(df.2, RUN), aes(x=variable, y=value), color="red", size=2)
dev.off()

# arrow plot

par(mfcol=c(3,4),mai=c(0.5,0.6,0.3,0.2))
plot(c(), xlim=c(0,1000000), ylim=c(0,1000000),xlab="NPOP1",ylab="NPOP2")
for (i in 2:(nrow(modelSPA))) {
  arrows(modelSPA$NPOP1[1], modelSPA$NPOP2[1], modelSPA$NPOP1[i], modelSPA$NPOP2[i], col="red")  }
plot(c(), xlim=c(0,1000000), ylim=c(0,1000000),xlab="NPOP3",ylab="NPO4")
for (i in 2:(nrow(modelSPA))) {
  arrows(modelSPA$NPOP3[1], modelSPA$NPOP4[1], modelSPA$NPOP3[i], modelSPA$NPOP4[i], col="red")  }
plot(c(), xlim=c(0,1000000), ylim=c(0,1000000),xlab="NANC1",ylab="NANC2")
for (i in 2:(nrow(modelSPA))) {
  arrows(modelSPA$NANC1[1], modelSPA$NANC2[1], modelSPA$NANC1[i], modelSPA$NANC2[i], col="red") }
plot(c(), xlim=c(0,1000000), ylim=c(0,1000000),xlab="NANC1",ylab="NANC3")
for (i in 2:(nrow(modelSPA))) {
  arrows(modelSPA$NANC1[1], modelSPA$NANC3[1], modelSPA$NANC1[i], modelSPA$NANC3[i], col="red")  }
plot(c(), xlim=c(100,10000), ylim=c(0,30000),xlab="TMRMG",ylab="TDIV1")
for (i in 2:(nrow(modelSPA))) {
  arrows(modelSPA$TMRMG[1], modelSPA$TDIV1[1], modelSPA$TMRMG[i], modelSPA$TDIV1[i], col="red") }
plot(c(), xlim=c(0,100000), ylim=c(0,100000),xlab="TDIV2",ylab="TDIV3")
for (i in 2:(nrow(modelSPA))) {
  arrows(modelSPA$TDIV2[1], modelSPA$TDIV3[1], modelSPA$TDIV2[i], modelSPA$TDIV3[i], col="red")  }
plot(c(), xlim=c(-0.0005,0.02), ylim=c(-0.0005,0.02),xlab="MIG01R",ylab="MIG10R")
for (i in 2:(nrow(modelSPA))) {
  arrows(modelSPA$MIG01R[1], modelSPA$MIG10R[1], modelSPA$MIG01R[i], modelSPA$MIG10R[i], col="red") }
plot(c(), xlim=c(-0.0005,0.02), ylim=c(-0.0005,0.02),xlab="MIG12R",ylab="MIG21R")
for (i in 2:(nrow(modelSPA))) {
  arrows(modelSPA$MIG12R[1], modelSPA$MIG21R[1], modelSPA$MIG12R[i], modelSPA$MIG21R[i], col="red")  }
plot(c(), xlim=c(-0.0005,0.02), ylim=c(-0.0005,0.02),xlab="MIG23R",ylab="MIG32R")
for (i in 2:(nrow(modelSPA))) {
  arrows(modelSPA$MIG23R[1], modelSPA$MIG32R[1], modelSPA$MIG23R[i], modelSPA$MIG32R[i], col="red")  }
plot(c(), xlim=c(-0.005,0.01), ylim=c(-0.005,0.01),xlab="MIG01R",ylab="GR2")
for (i in 2:(nrow(modelSPA))) {
  arrows(modelSPA$MIG01R[1], modelSPA$GR2[1], modelSPA$MIG01R[i], modelSPA$GR2[i], col="red") }
```

---
### Blue regions

To identify sites at are more than expected (blue in compSFS)
CEU=0.7-1 hooded=0.0-0.1 (diff=>0.6)
Spain=1, hooded=0 and Spain=0, hooded=1 (diff>0.9)
```
paste 134inds_overlapped_filtered_norepeats_hwe_neuall_incIRQ_cor2to3_cnx1to3_keepcor2to3.frq 134inds_overlapped_filtered_norepeats_hwe_neuall_incIRQ_cor2to3_cnx1to3_keepcnx1to3.frq | awk 'NR>1 {print $1, $2, $5, $6, $11, $12}' | awk -v FS=OFS="\t|:" '{print $1, $2,$3,$4, $5, $6, $7, $8, $9, $10}' | awk '{print $1, $2, $4, $6, $8, $10}' | awk 'BEGIN{print "scaffold pos pop1_f1 pop1_f2 pop2_f1 pop2_f2"}1' > ceu_hooded_sfs.frq

paste 134inds_overlapped_filtered_norepeats_hwe_neuall_incIRQ_cor1_cnx1to3_keepcor1.frq 134inds_overlapped_filtered_norepeats_hwe_neuall_incIRQ_cor1_cnx1to3_keepcnx1to3.frq | awk 'NR>1 {print $1, $2, $5, $6, $11, $12}' | awk -v FS=OFS="\t|:" '{print $1, $2,$3,$4, $5, $6, $7, $8, $9, $10}' | awk '{print $1, $2, $4, $6, $8, $10}' | awk 'BEGIN{print "scaffold pos pop1_f1 pop1_f2 pop2_f1 pop2_f2"}1' > spain_hooded_sfs.frq

awk 'NR>1 {if ($3-$5>0.6||$4-$6>0.6) print $1, $2, $3, $5, $3-$5}' ceu_hooded_sfs.frq | awk '{if($NF ~ /^-/ ){sub(/-/, "", $NF); print $0}}' | awk 'BEGIN{print "scaffold pos pop1_f1 pop2_f1 pop2_f1, diff_af"}1' > ceu_hooded_sfs_above0.6.frq

awk 'NR>1 {if ($3-$5>0.9||$4-$6>0.9) print $1, $2, $3, $5, $3-$5}' spain_hooded_sfs.frq | awk '{if($NF ~ /^-/ ){sub(/-/, "", $NF); print $0}}' | awk 'BEGIN{print "scaffold pos pop1_f1 pop2_f1 pop2_f1, diff_af"}1' > spain_hooded_sfs_above0.9.frq

wc -l ceu_hooded_sfs_above0.6.frq #293
wc -l spain_hooded_sfs_above0.9.frq #1117

sort -k5 -r ceu_hooded_sfs_above0.6.frq 
```
The below scaffolds are found in the blue regions of ceu_hooded & spain_hooded, higher freq in hooded for both
scaffold_29: chr15 (73/293) ~335700bp 919873-1255588
scaffold_7: chr8 (161/293) ~150000bp 50858-198375

#### blast
```
ref="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2/genome_HC_allpaths41687_v2.5_chrW.fasta"

cd /dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2/fastsimcoal2/4Pop/cor1_cor2to3_cnx1to3_cnx6_50ind_folded_all_nochr18/fastsimcoal2/bestruns/blueregion

bedtools getfasta -fi /dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2/genome_HC_allpaths41687_v2.5_chrW.fasta -bed scaffold29.bed -fo scaffold29.fa

bedtools getfasta -fi /dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2/genome_HC_allpaths41687_v2.5_chrW.fasta -bed scaffold7.bed -fo scaffold7.fa
```
---

# Twisst
Tree-based approach to identify region(s) of introgression using various python scripts developed by Simon Martin (clone this): https://github.com/simonhmartin/genomics_general 

Added `3pop.pop` and `4pop.pop` for population-specific VCF. Added `3pop_groups_edited.tsv` and `4pop_groups_edited.tsv` for twisst step. <i> C. moneduloides</i> is used as an outgroup here. 

One tree generated for each block of 50 SNPs. 
No. of trees: 285,323 trees w/o scaffold 78+60, 287,217 with scaff78+60

Run `1.5.7_twisst.sh`

```
#!/bin/bash -l
#SBATCH -J twisst
#SBATCH --get-user-env
#SBATCH --mail-user=gwee@biologie.uni-muenchen.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=6
#SBATCH --time=2-00:00:00
#SBATCH -o /dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2/slurms/slurm-%j-%x.out

# sbatch 1.5.7_twisst.sh all 134inds_overlapped_filtered_norepeats_hwe_outgroup_biallele 4pop "-g cnx6 -g cnx3 -g cor1 -g cor2 -g O --outgroup O"
# sbatch 1.5.7_twisst.sh all 134inds_overlapped_filtered_norepeats_hwe_outgroup_biallele 3pop "-g cnx3 -g cor1 -g cor2 -g O --outgroup O"

module load vcftools
module load bcftools
conda activate biotools
#picard.jar=2.25.7
#module load plink2
#using plink1.9 instead

ref="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2/genome_HC_allpaths41687_v2.5_chrW.fasta"
dat="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2"
scaff="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/01_probes/scaffold2chr"
gentools="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2/06_results/fst/genomics_general"

echo $(date)
STARTTIME=$(date +%s)

cd $dat/06_results/twisst

#subsample for twisst
vcftools --keep $3.pop --gzvcf ${dat}/05.1_recal/overlap/${2}.vcf.gz --recode --out $1_$3
mv $1_$3.recode.vcf $1_$3.vcf
bgzip $1_$3.vcf
bcftools index $1_$3.vcf.gz
bcftools view --threads 6 -V indels --min-alleles 2 -Oz -o $1_$3_filtered.vcf.gz $1_$3.vcf.gz
bcftools view -r $(echo $(cat chr.txt) | sed 's/ //g') -Oz -o $1_$3_filtered_100scaff.vcf.gz $1_$3.vcf.gz

#Phase with Beagle5.4
java -Xmx12g -jar beagle.jar gt=$1_$3_filtered_100scaff.vcf.gz out=$1_$3_filtered_100scaff_phased impute=true nthreads=6 window=1 overlap=0.1
python ${gentools}/VCF_processing/parseVCF.py -i $1_$3_filtered_100scaff_phased.vcf.gz | gzip > $1_$3_filtered_phased.geno.gz

#NJ tree per 50 SNPs
python ${gentools}/phyml_sliding_windows.py -T 6 -g $1_$3_filtered_phased.geno.gz --prefix $1_$3.phyml_bionj.w50 -w 50 --windType sites --model GTR --optimise n

awk '{print $1"_A", $2}' $3_groups.tsv > $3_groups_editedA.tsv
awk '{print $1"_B", $2}' $3_groups.tsv > $3_groups_editedB.tsv
cat $3_groups_editedA.tsv $3_groups_editedB.tsv | awk '{print $2, $1}' | sort | awk '{print $2, $1}' OFS="\t" > $3_groups_edited.tsv
rm $3_groups_editedA.tsv $3_groups_editedB.tsv

python twisst.py -t $1_$3.phyml_bionj.w50.trees.gz -w $1_$3.weights.csv.gz --method complete --groupsFile $3_groups_edited.tsv $4
zcat $1_$3.weights.csv.gz | sed "s/outgroup/O/g" >  $1_$3.weights.edited.csv
gzip $1_$3.weights.edited.csv
```
---

# Fastsimcoal to Twisst

Simulating the best-fitted estimates for model7x and model8x each to trees and run with Fastsimcoal. 

See folder `fastsimcoal/4Pop/cor1_cor2to3_cnx1to3_cnx6_50ind_folded_all_nochr18/fastsimcoal2/twisst` for input files

Edit the maxL.par file of each model. Remove IRQ and replace with an outgroup C. moneduloides to prepare 3pop (excluding outgroup) trees. Assuming the outgroup diverged 10M years ago from C. corone, the generation time is 1727115. Specify to simulate 285323 trees of 2500 bp each. A total of 14021570 and 14577480 polymorphic sites were generated from 7x and 8x, respectively (mean no. of SNPs per tree = 49 and 51, respectively). 

```
//Number of independent loci [chromosome] 
285323 0
//Per chromosome: Number of contiguous linkage Block: a block is a set of contiguous loci
1
//per Block:data type, number of loci, per generation recombination and mutation rates and optional parameters
DNA 2500 0 3.18e-9 0.5
```

Generate trees with `1.6.6_fastsimcoalTwisst.sh`
```
# sbatch 1.6.6_fastsimcoalTwisst.sh 3PopModel7x 3popsim run1
# sbatch 1.6.6_fastsimcoalTwisst.sh 3PopModel8x 3popsim run1

module load vcftools/0.1.14-gcc8
module load bcftools/1.10.2-gcc8
conda activate easySFS

ref="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2/genome_HC_allpaths41687_v2.5_chrW.fasta"
dat="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2"

echo $(date)
STARTTIME=$(date +%s)

cd ${dat}/fastsimcoal2/4Pop/cor1_cor2to3_cnx1to3_cnx6_50ind_folded_all_nochr18/fastsimcoal2/bestruns/twisst

#produce treefile
cp ${1}tree.par ${1}tree_${3}.par
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
```
## Plotting for Twisst results
To plot empirical trees with and without barrier locus and also simulated trees 
with R and ggtern

`plot_all_3pop.R` & `plot_all_3pop_simulationcomparison.R`
```
library(ape, lib.loc="/dss/dsshome1/lxc0E/di67kah/R")
library(ggtern)
library(ggplot2)
library(dplyr)

setwd("/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2/06_results/twisst")
source("plot_twisst.R")

hash = paste("#topo1 (O,((EURn,SPA),EURc));",
             "#topo2 (O,((EURn,EURc),SPA));",
             "#topo3 (O,(EURn,(EURc,SPA)));", sep="\n")

weights_file <- "all_3pop.weights.edited.csv.gz"
window_data_file <- "all_3pop.phyml_bionj.w50.data.tsv"
w1 <- read.csv(weights_file, skip = 3, header=TRUE, sep = "\t")
t1 <- read.csv(window_data_file, header=TRUE, sep = "\t")
dat <- cbind(t1,w1)


writeLines(hash, con = "all_3pop.weights.edited.name.csv")
write.table(dat[7:length(dat)],file="all_3pop.weights.edited.name.csv", sep= "\t", quote=FALSE, append=TRUE)
system("gzip all_3pop.weights.edited.name.csv")

twisst_data <- import.twisst(weights_files="all_3pop.weights.edited.name.csv.gz", window_data_files=window_data_file)
pdf("all_3pop_plotsummary_edited.pdf") #6x14
plot.twisst.summary(twisst_data, lwd=3, cex=1)
dev.off()

######################### identifying scaffold78 and 60 region ###############################

dat1 <- dat %>% group_by(scaffold) %>%
  filter(scaffold != "scaffold_78" & scaffold != "scaffold_60")

w1edit <- dat1[7:length(dat)]
t1edit <- dat1[1:6]

writeLines(hash, con = "all_3pop.weights.edited.nochr18.csv")
write.table(w1edit,file="all_3pop.weights.edited.nochr18.csv", sep= "\t", quote=FALSE, append=TRUE)
system("gzip all_3pop.weights.edited.nochr18.csv")
write.table(t1edit,file="all_3pop.phyml_bionj.w50.data.nochr18.tsv", sep= "\t", quote=FALSE)

weights_file_nochr18 <- "all_3pop.weights.edited.nochr18.csv.gz"
window_data_file_nochr18 <- "all_3pop.phyml_bionj.w50.data.nochr18.tsv"
twisst_data_nochr18 <- import.twisst(weights_files=weights_file_nochr18, window_data_files=window_data_file_nochr18)
plot.twisst.summary(twisst_data_nochr18, lwd=3, cex=1)

dat2 <- dat %>% group_by(scaffold) %>%
  filter(scaffold == "scaffold_78" | scaffold == "scaffold_60")

w2edit <- dat2[7:length(dat)]
t2edit <- dat2[1:6]

writeLines(hash, con = "all_3pop.weights.edited.chr18.csv")
write.table(w2edit,file="all_3pop.weights.edited.chr18.csv", sep= "\t", quote=FALSE, append=TRUE)
system("gzip all_3pop.weights.edited.chr18.csv")
write.table(t2edit,file="all_3pop.phyml_bionj.w50.data.chr18.tsv", sep= "\t", quote=FALSE)

weights_file_chr18 <- "all_3pop.weights.edited.chr18.csv.gz"
window_data_file_chr18 <- "all_3pop.phyml_bionj.w50.data.chr18.tsv"
twisst_data_chr18 <- import.twisst(weights_files=weights_file_chr18, window_data_files=window_data_file_chr18)
plot.twisst.summary(twisst_data_chr18, lwd=3, cex=1)


########################### ternary plot ###########################
dattern <- mutate(dat, type = case_when (scaffold == "scaffold_78" | scaffold == "scaffold_60" ~ "barrier",
                                         scaffold != "scaffold_78" | scaffold != "sacffold_60" ~ "neutral"))

#### MAIN ####
pdf("ternaryplot_3pop_all.pdf")
ggtern(data=dattern,aes(x=topo1,y=topo2, z=topo3))+ geom_hex_tern(bins=150) + 
  scale_shape_manual(values=c(15, 17)) +
  scale_fill_distiller(palette = "RdPu") +
  #scale_fill_gradient(low = "pink", high = "darkred")  +
  scale_L_continuous(breaks = 0:5 / 5, labels = 0:5/ 5) +
  scale_R_continuous(breaks = 0:5 / 5, labels = 0:5/ 5) +
  scale_T_continuous(breaks = 0:5 / 5, labels = 0:5/ 5) +
  guides(fill = guide_colorbar(order = 1), alpha = guide_none()) +
  theme_custom(base_size = 14, base_family = "", tern.plot.background = NULL,
               tern.panel.background = "white", col.T = "limegreen", col.L = "royalblue", col.R = "darkorange", col.grid.minor = "white") +
  theme_showarrows() + 
  #theme_noarrows() +
  theme_legend_position('topleft') +
  labs(fill   = "Density of sites",
       Tarrow = "ILS + Introgression",
       Larrow = "ILS",
       Rarrow = "Barrier")
dev.off()
##only neutral sites
emp7 <- ggtern(dattern %>% filter(type == "neutral"),aes(x=topo1,y=topo2, z=topo3))+ geom_hex_tern(bins=150) + 
  scale_shape_manual(values=c(15, 17)) +
  scale_fill_distiller(palette = "RdPu") +
  #scale_fill_gradient(low = "pink", high = "darkred")  +
  scale_L_continuous(breaks = 0:5 / 5, labels = 0:5/ 5) +
  scale_R_continuous(breaks = 0:5 / 5, labels = 0:5/ 5) +
  scale_T_continuous(breaks = 0:5 / 5, labels = 0:5/ 5) +
  guides(fill = guide_colorbar(order = 1), alpha = guide_none()) +
  theme_custom(base_size = 14, base_family = "", tern.plot.background = NULL,
               tern.panel.background = "white", col.T = "limegreen", col.L = "royalblue", col.R = "darkorange", col.grid.minor = "white") +
  theme_showarrows() + 
  #theme_noarrows() +
  theme_legend_position('topleft') +
  labs(fill   = "Density of sites",
       Tarrow = "ILS + Introgression",
       Larrow = "ILS",
       Rarrow = "Barrier") +
  ggtitle("Observed neutral sites") +
  theme(plot.title = element_text(hjust = 0.5, size = 15))
pdf("ternaryplot_3pop_neutral.pdf")
emp7
dev.off()

#### ternary contour ####
## conducting density estimation on barrier (scaffold78+60) and neutral sites separately
c1 <- ggtern(dattern, aes(x = topo1,y = topo2, z = topo3, color = type)) +
  stat_density_tern(aes(alpha = ..level.., fill = type), 
                    geom = 'polygon', 
                    bins = 50,
                    color = "grey",bdl=0.005) +
  #geom_point(alpha = 0.4, size=1) +
  scale_fill_manual(values = c("#DC3220", "#005AB5")) +
  theme_custom(base_size = 14, base_family = "", tern.plot.background = NULL,
               tern.panel.background = "white", col.T = "limegreen", col.L = "royalblue", col.R = "darkorange", col.grid.minor = "white") +
  theme_showarrows() + 
  #theme_noarrows() +
  theme_legend_position('topleft') +
  labs(fill   = "Type of sites",
       alpha = "Density",
       Tarrow = "ILS + Introgression",
       Larrow = "ILS",
       Rarrow = "Barrier") +
  theme_hidegrid_major()
## conducting density estimation on all loci collectively
c2 <- ggtern(dattern, aes(x = topo1,y = topo2, z = topo3)) +
  stat_density_tern(aes(alpha = ..level..), 
                    geom = 'polygon', 
                    bins = 50,
                    color = "grey",bdl=0.005) +
  scale_fill_manual(values = c("#DC3220", "#005AB5")) +
  theme_custom(base_size = 14, base_family = "", tern.plot.background = NULL,
               tern.panel.background = "white", col.T = "limegreen", col.L = "royalblue", col.R = "darkorange", col.grid.minor = "white") +
  theme_showarrows() + 
  theme_legend_position('topleft') +
  labs(fill = "Type of sites",
       alpha = "Density",
       Tarrow = "ILS + Introgression",
       Larrow = "ILS",
       Rarrow = "Barrier") +
  theme_hidegrid_major()

pdf("ternarycontour_3pop.pdf",width = 14, height = 8.5)
grid.arrange(c2,c1,ncol=2)
dev.off()

############ simulated twisst ############
setwd("/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2/fastsimcoal2/4Pop/cor1_cor2to3_cnx1to3_cnx6
_50ind_folded_all_nochr18/fastsimcoal2/bestruns/twisst")
weights_file_3popsim7 <- "./3PopModel7xtree_run2/3PopModel7x_run2.weights.csv.gz"
weights_file_3popsim8 <- "./3PopModel8xtree_run2/3PopModel8x_run2.weights.csv.gz"
window_data_file_3popsim7 <- "./3PopModel7xtree_run2/3PopModel7x_run2.win.data.tsv"
window_data_file_3popsim8 <- "./3PopModel8xtree_run2/3PopModel8x_run2.win.data.tsv"

#edit files due to name changes
w3popsim7 <- read.csv(weights_file_3popsim7, skip = 3, header=TRUE, sep = "\t")
colnames(w3popsim7) <- c("topo3","topo1","topo2")
w3popsim7 <- relocate(w3popsim7,topo3, .after = topo2)
t3popsim7 <- read.csv(window_data_file_3popsim7, header=TRUE, sep = "\t")
dat3popsim7 <- cbind(t3popsim7 ,w3popsim7)
w3popsim8 <- read.csv(weights_file_3popsim8, skip = 3, header=TRUE, sep = "\t")
colnames(w3popsim8) <- c("topo3","topo1","topo2")
w3popsim8 <- relocate(w3popsim8,topo3, .after = topo2)
t3popsim8 <- read.csv(window_data_file_3popsim8, header=TRUE, sep = "\t")
dat3popsim8 <- cbind(t3popsim8 ,w3popsim8)

sim7 <- ggtern(dat3popsim7,aes(x=topo1,y=topo2, z=topo3))+ geom_hex_tern(bins=150) + 
  scale_shape_manual(values=c(15, 17)) +
  scale_fill_distiller(palette = "RdPu") +
  #scale_fill_gradient(low = "pink", high = "darkred")  +
  scale_L_continuous(breaks = 0:5 / 5, labels = 0:5/ 5) +
  scale_R_continuous(breaks = 0:5 / 5, labels = 0:5/ 5) +
  scale_T_continuous(breaks = 0:5 / 5, labels = 0:5/ 5) +
  guides(fill = guide_colorbar(order = 1), alpha = guide_none()) +
  theme_custom(base_size = 14, base_family = "", tern.plot.background = NULL,
               tern.panel.background = "white", col.T = "limegreen", col.L = "royalblue", col.R = "darkorange", col.grid.minor = "white") +
  theme_showarrows() + 
  #theme_noarrows() +
  theme_legend_position('topleft') +
  labs(fill = "Density",
    Tarrow = "ILS + Introgression",
    Larrow = "ILS",
    Rarrow = "Barrier")+
  ggtitle(paste("Simulated neutral sites from", "genome-wide introgression model", sep = "\n"))+
  theme(plot.title = element_text(hjust = 0.5, size = 15))

sim8 <- ggtern(dat3popsim8,aes(x=topo1,y=topo2, z=topo3))+ geom_hex_tern(bins=150) + 
  scale_shape_manual(values=c(15, 17)) +
  scale_fill_distiller(palette = "RdPu") +
  #scale_fill_gradient(low = "pink", high = "darkred")  +
  scale_L_continuous(breaks = 0:5 / 5, labels = 0:5/ 5) +
  scale_R_continuous(breaks = 0:5 / 5, labels = 0:5/ 5) +
  scale_T_continuous(breaks = 0:5 / 5, labels = 0:5/ 5) +
  guides(fill = guide_colorbar(order = 1), alpha = guide_none()) +
  theme_custom(base_size = 14, base_family = "", tern.plot.background = NULL,
               tern.panel.background = "white", col.T = "limegreen", col.L = "royalblue", col.R = "darkorange", col.grid.minor = "white") +
  theme_showarrows() + 
  #theme_noarrows() +
  theme_legend_position('topleft') +
  labs(fill   = "Density",
       Tarrow = "ILS + Introgression",
       Larrow = "ILS",
       Rarrow = "Barrier")+
  ggtitle(paste("Simulated neutral sites from", "locus-specific introgression model",sep="\n")) +
  theme(plot.title = element_text(hjust = 0.5, size = 15))
pdf("ternaryplot_emp_sim7_sim8.pdf",width=16,height=6)
grid.arrange(emp7,sim7,sim8, nrow=1)
dev.off()

#### ternary count for neutral sites ####
p1 <- ggtern(dattern %>% filter(type == "neutral"),aes(topo1,topo2,topo3)) + 
  theme_bw() +
  geom_tri_tern(bins=20,aes(fill=after_stat(count))) + 
  stat_tri_tern(bins=20,
                geom='text',
                aes(label=sprintf("%.0f",after_stat(count))),
                size=2, color='white',centroid=TRUE)
p2 <- ggtern(dat3popsim7,aes(topo1,topo2,topo3)) + 
  theme_bw() +
  geom_tri_tern(bins=20,aes(fill=after_stat(count))) + 
  stat_tri_tern(bins=20,
                geom='text',
                aes(label=sprintf("%.0f",after_stat(count))),
                size=2, color='white',centroid=TRUE)
p3 <- ggtern(dat3popsim8,aes(topo1,topo2,topo3)) + 
  theme_bw() +
  geom_tri_tern(bins=20,aes(fill=after_stat(count))) + 
  stat_tri_tern(bins=20,
                geom='text',
                aes(label=sprintf("%.0f",after_stat(count))),
                size=2, color='white',centroid=TRUE)

## using the count to compare against simulated results
emp <- ggplot_build(p1)$data[[1]] %>% group_by(group)  %>% 
  filter(row_number()==1) %>%
  dplyr::mutate(count = replace_na(count, 0))
sim7 <- ggplot_build(p2)$data[[1]] %>% group_by(group)  %>% 
  filter(row_number()==1) %>%
  dplyr::mutate(count = replace_na(count, 0))
sim8 <- ggplot_build(p3)$data[[1]] %>% group_by(group)  %>% 
  filter(row_number()==1) %>%
  dplyr::mutate(count = replace_na(count, 0))

comp <- cbind(emp, "sim7"=sim7$count, "sim8"=sim8$count) %>%
  mutate(sim7diff=(sim7-count), sim8diff=(sim8-count))

sim7diff <- data.frame(comp$x,comp$y,comp$z,comp$sim7diff) %>%
  #uncount(abs(comp$sim7diff)) %>%
  rename("x"="comp.x","y"="comp.y","z"="comp.z","count"="comp.sim7diff")
sim8diff <- data.frame(comp$x,comp$y,comp$z,comp$sim8diff) %>%
  rename("x"="comp.x","y"="comp.y","z"="comp.z","count"="comp.sim8diff")

p7 <- ggtern(data = sim7diff, aes(x=x,y=y,z=z)) +
  geom_point(size=4, aes(color=count)) +
  theme_bw() +
  scale_colour_gradient2(low = muted("darkred"),mid = "white",high = muted("darkblue"),midpoint = 0, limits = c(-10000, 10000)) +
  theme_custom(base_size = 14, base_family = "", tern.plot.background = NULL,
               tern.panel.background = "white", col.T = "limegreen", col.L = "royalblue", col.R = "darkorange", col.grid.minor = "white") +
  theme_showarrows() + 
  theme_legend_position('topleft') +
  labs(color  = "Count",
       Tarrow = "ILS + Introgression",
       Larrow = "ILS",
       Rarrow = "Barrier")  +
  ggtitle("Exp(GWS)-Obs") +
  theme(plot.title = element_text(hjust = 0.5, size = 15))

p8 <- ggtern(data = sim8diff, aes(x=x,y=y,z=z)) +
  geom_point(size=4, aes(color=count)) +
  theme_bw() +
  scale_colour_gradient2(low = muted("darkred"),mid = "white",high = muted("darkblue"),midpoint = 0, limits = c(-10000, 10000)) +
  theme_custom(base_size = 14, base_family = "", tern.plot.background = NULL,
               tern.panel.background = "white", col.T = "limegreen", col.L = "royalblue", col.R = "darkorange", col.grid.minor = "white") +
  theme_showarrows() + 
  theme_legend_position('topleft') +
  labs(fill   = "Type of sites",
       alpha = "Density",
       Tarrow = "ILS + Introgression",
       Larrow = "ILS",
       Rarrow = "Barrier")  +
  ggtitle("Exp(LSI)-Obs") +
  theme(plot.title = element_text(hjust = 0.5, size = 15))
```

---

# Jaatha
written by Dirk Metzler

In addition to fastsimcoal2 we used the composite-likelihood optimization
method jaatha (Mathew et al. 2013, https://doi.org/10.1002/ece3.722)
to optimize to estimate demographic model parameters. 
This allowed us to complement folded joint site-frequency spectra with additional
statistics, aiming to improve estimations of demographic parameters and in
addition to estimate the recombination rate.

## Parameter ranges

The options for parameter scaling differ between fastsimcoal2 and jaatha.
For the two methods we were able to choose parameter ranges that vastly overlap
and both contain the plausible ranges but differ a bit at their boundaries. 
The parameter ranges of the models used in jaatha were as follows:

| parameter                                  |   lower limit     |     upper limit      |
| ------------------------------------------ | ----------------- | -------------------- |
|any present population size                 |     100           |     500000           | 
|any ancestral population size               |     100           |     500000           |
|generations between population split events |     100           |     100000           |
|growth rate per generation in $W$ or $H$    |       0           |     0.0001           |
|migration rates per generation              |       0           |     0.001            |
|generations since migration between $W$ and $H$ started  |  100 |     10000            |
|recombination rate                          |      $\mu/3$      |   $30\cdot\mu$       |

Here, $\mu$ is the mutation rate, for which we assumed the value $3.18\cot10^{-9}$, following __ADD REFERENCE__.
Population growth in $W$ and $H$ was allowed to start at any time point since the most
recent split of the repspective population from a joint ancestral population with its
sister population.

## Statistics used in the jaatha analysis

For each pair of populations $S$—$W$, $S$—$H$, $S$—$I$, $W$—$H$, $W$—$I$, and $H$—$I$
we used the coarsening of the joint site-site frequency spectrum as described in
Naduvilezhath, Rose, Metzler (2011, https://doi.org/10.1111/j.1365-294X.2011.05131.x),
see also Tellier et al. (2011, https://doi.org/10.1371/journal.pone.0018155)
resulting in 14 categories of polymorphism for each of the six population pairs.

In addition to the frequency-spectrum based statisitcs
we used statistics that are based on violations of the four-gamete
condition in pairs of polymorphic sites (Hudson, Caplan,
https://doi.org/10.1093/genetics/111.1.147).
The four-gamete condition is violated for two polymorphic sites if there is no
tree that explains both polymorphic sites with only one mutation per site.
This may occur due to a recombination event between the sites or if one of the
sites was hit by more than one mutation.
The four-gamete condition is violated if and only if haplotypes of all four
combinations of alleles exist.
As we apply this to unphased data, we count the four-gamete condition as
violated only if it is violated in every possible phasing of the data from the two
sites.
For 91 genomic ranges of length 10000 bp sequenced for 15 crows from $S$, 24 crows from $W$
54 crows from $H$ and 5 crows from $I$
we calculated the fraction of pairs of sites that violated the four-gamete condition among
all pairs of sites with polymorphism for which each allele existed in two or more sequences.
Further for the pairs populations $S$—$W$, $W$—$H$ and $H$—$I$ we calculated the
fraction of pairs or sites that violated the four gamete condition only if sequences of
both populations were pooled, again relative to the number of pairs of sites in which
both alleles existed in two or more sequences.
For each such fraction we categorized the 91 genomic regions into four groups.
For the single populations the grouping criterion was whether the fraction was in
[0, 0.1), [0.1, 0.2), [0.2, 0.3] or [0.3, 1], and for the population pairs the
intervals were [0, 0.005), [0.005, 0.01), [0.01, 0.015) and [0.015, 1].
The statistics used in addition to frequency-spectra based statistics
in our analysis are then the 28 numbers of genomic ranges falling into
each of the four categories for the seven fractions from four single populations and three
pairs of populations.

## Simulation and optimization strategies

To simulate frequency spectra for parameter combinations during the jaatha runs, we used msprime
(Baumdicker et al., 2022, https://doi.org/10.1093/genetics/iyab229),
simulating 400 loci of 2000 bp without recombination at a time.
The observed numbers were then scaled up to approximate the expeced frequency spectra for full genome data.
To simulate 91 loci of 10000 bp with recombination to approximate the expected values for
the four-gamete condition based statistics, we used the R package scrm
(Staab et al. 2015, https://doi.org/10.1093/bioinformatics/btu861).
In each iteration step of jaatha, 1000 such simulations were performed.
To simulate mutations within jaatha we combined both msprime and scrm with the infinite sites model.
For the search of the initial parameter values in jaatha we applied the zoom-in method and then
repeated the search for an optimum five times.

## Bootstrap simulations

To assess the accuracy of our fastsimcoal2 and jaatha estimations, we applied a parametric bootstrap
approach. For this, we launched 100 simulations with using the fitted parameters.
We repeated this in different variants, both with infinite-sites simulations and with finite-sites simulations
to explore a possible bias by neglecting back mutations and double hits.
For the bootstrap simulations we simulated data of a similar size as the original ones in msprime,
that is 912054 loci of 1000 bp.
For the finite-sites bootstrap simulations we explored two different variants,
one combining msprime and scrm with seq-gen
(Rambaut, A. and Grassly, N. C. (1997) Seq-Gen: An application for the Monte Carlo simulation of
DNA sequence evolution along phylogenetic trees. Comput. Appl. Biosci. 13: 235-238.)
with the HKY model with a transition–transversion ratio of 2
(Hasegawa, Kishino, Yano, 1985, https://doi.org/10.1007/BF02101694),
and the other using the finite-sites simulator for the Jukes-Cantor model (__TRUE?__)
implemented in fastsimcoal2.
Also for infinite-sites simulations we applied both the msprime/scrm simulations and
the infinite-sites simulator avaialable in fastsimcoal.
The re-estimates were then carried out with the same procedures as for the orignal data,
that is with infinite-sites mutations
and upscaling in of the frequency-spectrum simulations in the case of jaatha.
For any parameter that was origninally estimated as $\hat\theta$, and for which $\bar\theta$ was the mean
of the bootstrap estimations and $q_{0.025}$ and $q_{0.975}$ are their 0.025 and 0.975 quantiles,
we consider $\hat\theta-(\bar\theta-\hat\theta)=2\hat\theta-\bar\theta$ as bootstrap bias-correced
estimation and $[2\hat\theta-q_{0.975},\ 2\hat\theta-q_{0.025}]$ as boostrap confidence interval
(Efron, Tibshirani; An Introduction to the Bootstrap). Quantiles were calculated with default settings
of the quantile funcition in R.

## Results

### $W$ split off from $S$ or $H$?

In the best fitting parametrizations of our demographic models,
population $W$, consisting of the black carrion crows outside of the
Iberian Peninsula, first split off from the Iberian population $S$ before
gene flow between $W$ and the hooded crow population $H$ set in.
In contrast to this, the parameter estimations obtained with the restriction
that $W$ initially split off from $H$ obtained a substantially smaller
composite log-likelihood (.... vs. ... for the fastsimcoal2 estimations
and  -5078795 vs. -4350459 for the jaatha estimations) and predictied a
substantially higher number of fixed differences between $S$ and $W$ and
polymorphisms in which one allele has a high frequency in $S$ and the other
allele has a high frequency in $W$, see Figures...., red areas in top left
and bottom right of plot...
Even though the estimations obtained with fastsimcoal2 and jaatha differ
slightly in the precise parameter values, they agree in the above points.

### Migration rates into $W$ and $S$—$W$ split time.

According to the jaatha estimation, $W$ split off from the $S$ population
6270 generations ago and after that, 0.0117 % of the crows of each generation 
in $W$ were immigrants from $S$. Further, 0.0287 % of the
crows in $W$ in each of the last 5924 generations were immigrants from $H$.
This implies that 90.9 % of the ancestral lineages in $W$ trace back their
ancestry into the $S$ or the $H$ population with immigration during this time span
and 70.9 % of them, that is 64.5 % of all lineages, come from the hooded crow population $H$
(not taking into account the past of lineages before their most recent
migration events).

### Growth rates

According to the jaatha estimation, $W$ grew during the last 3500
generation with a rate of $4.9\cdot 10^{-5}$ per generation, thus increasing
in size by 18.9 %. Accoring to the same model fit, the hooded crows
population $H$ began growing already approx. 12650 generation ago and has
increased in size by a factor of 2.35 since then.

### Accuracy of demographic parameter inferrence

### Inspection of site frequency spectra

As models in general, our demographic model is still a simplification of
reality and for the joint site-frequency spectrum of $W$ and $H$, for example,
predicts more polymorphisms in which one allele has a frequency below 0.2
in $W$ and between 0.6 and 0.8 in $H$ than empirically observered and in turn
fewer in which one allele has a frequency between 0.6 and 0.8 in $W$ and below
0.2 in $H$, see figures/compfrspemp.pdf red areas in the top left and
blue area in the bottom right of plot WH.
[For the caption of this and similar plots: if $a$ is the observered and $b$ the
simulated number of polymorphisms of certain frequencies, the plots show $a/(a+b)$
on a scale from red to blue; red means there were observed but no simulated
polymorphisms, blue means there were simulated but no observed polymorphism and
white means that the numbers of observed and predicted polymorphisms were equal.]
(For better comparison with the first case, note that the latter is equivalent
to one allele having a frequency between 0.2 and 0.4 in $W$ and above 0.8 in
$H$.)
This was similar when infinite-sites simulation were used to generate frequency
spectra for comparison (figures/compfrspemp_is.pdf) and also for simulations
(with finite sites) with best parameters found for the model in which $W$ first
split off from the $H$ population (figure/compfrhdemp.pdf).
In the latter case, the worse fit of that model becomes apparent from the red
regions in plot SW as the model would predict more polymorphisms for which one
allele is fixed or at least of high frequency in $S$ and the other is fixed or
of high freqency in $W$.
The finite-sites simulation for model 7x fitted with fastsimcoal2 showed
a similar pattern as the simulations with the best models fitted with jaatha,
but with the differences for allels that are frequent in $W$ and rare in $H$ a
bit more pronounced (figure/compm7xemp.pdf) and generally a slight increase of
shared polymorphims of intermediate frequencies in both populations.

### Parametric bootstrapping

The parametric bootstrap analyses of our parameter estimations indicate for
all model parameters a certain amount of bias and uncertainty.
This applies especially to ancestral population sizes, which are hard to
estimate given high rates of recent gene flow between the present populations.
There is a slight difference between the bias of estimations for infinite-site
simulations and for finite-site simulations, but this inaccuracy adds only
little to the general uncertainty of the precise parameter estimations
(Figure figures/frsp_bs_is_fs.pdf; blue refers to boostrap results with
finite-sites simulations and black to those with infinite-sites simulations;
note that in Figures I and J the axes do not beginn at 0.)
The fit of pairwise site-frequency spectra simulated with bootstrap-based
bias-corrected parameter values to observed frequency spectra was similar to
that of parameter estimations without this correction
(figures/compfrspbsemp.pdf vs. figures/compfrspemp.pdf).

Table for supplement/appendix of parameter estimations and bootstrap corrected values and confidence ranges
(in msprime scaling):

|                     |     orig  |   bscor  |   upper  |   lower  | 
|---------------------|-----------|----------|----------|----------| 
| Ne S                |  8.33e+03 | 9.70e+03 | 1.11e+04 | 8.68e+03 |
| Ne W                |  3.41e+04 | 4.06e+04 | 4.86e+04 | 3.47e+04 |  
| Ne H                |  4.98e+04 | 5.67e+04 | 6.62e+04 | 4.93e+04 | 
| Ne I                |  1.25e+04 | 1.34e+04 | 1.48e+04 | 1.22e+04 | 
| Ne Anc S-W          |  4.70e+05 | 4.95e+05 | 6.29e+05 | 4.70e+05 | 
| Ne Anc H-I          |  4.58e+05 | 4.70e+05 | 5.69e+05 | 4.46e+05 | 
| Ne Anc S-W-H-I      |  2.96e+04 | 2.97e+04 | 3.03e+04 | 2.93e+04 | 
| t S-W               |  6.27e+03 | 7.67e+03 | 9.70e+03 | 6.27e+03 | 
| t H-I               |  1.34e+04 | 1.55e+04 | 2.13e+04 | 1.22e+04 | 
| t (S-W)-(H-I)       |  5.44e+04 | 5.46e+04 | 6.24e+04 | 4.95e+04 | 
| mig W <- S          |  1.17e-04 | 9.55e-05 | 1.12e-04 | 7.03e-05 | 
| mig S <- W          |  9.50e-05 | 7.58e-05 | 8.98e-05 | 5.93e-05 | 
| mig W <- H          |  2.87e-04 | 2.33e-04 | 2.76e-04 | 1.71e-04 | 
| mig H <- W          |  8.78e-05 | 7.29e-05 | 8.38e-05 | 6.06e-05 | 
| mig I <- H          |  6.63e-05 | 5.84e-05 | 6.55e-05 | 5.13e-05 | 
| mig H <- I          |  8.82e-06 | 8.82e-06 | 8.83e-06 | 8.76e-06 | 
| t start mig W-H     |  5.92e+03 | 7.14e+03 | 8.58e+03 | 6.05e+03 | 
| Growth rate Ne W    |  4.94e-05 | 6.78e-05 | 9.71e-05 | 1.99e-05 | 
| Start time growth W |  3.50e+03 | 3.87e+03 | 9.96e+03 | 8.63e+02 | 
| Growth rate Ne H    |  6.75e-05 | 8.07e-05 | 1.29e-04 | 3.75e-05 | 
| Start time growth H |  1.27e+04 | 2.08e+04 | 3.66e+04 | 1.17e+04 | 
| recomb rate         |  3.20e-09 | 3.21e-09 | 3.23e-09 | 3.07e-09 | 

same with fastsimcoal scaling (sign or growth rate is different) and parameter naming:

|         |      orig  |     bscor  |     lower  |     upper | 
|---------|------------|------------|------------|-----------| 
| NANC3   | 5.924e+04  | 5.948e+04  | 5.851e+04  | 6.059e+04 | 
| TDIV3   | 5.442e+04  | 5.463e+04  | 4.948e+04  | 6.241e+04 | 
| NANC1   | 9.408e+05  | 9.894e+05  | 9.402e+05  | 1.258e+06 | 
| NANC2   | 9.168e+05  | 9.392e+05  | 8.928e+05  | 1.138e+06 | 
| TDIV2   | 1.338e+04  | 1.550e+04  | 1.217e+04  | 2.133e+04 | 
| NPOP3   | 9.967e+04  | 1.135e+05  | 9.854e+04  | 1.324e+05 | 
| GROW2   | -6.748e-05 | -8.067e-05 | -3.746e-05 | -1.293e-04 | 
| NPOP4   | 2.509e+04  | 2.680e+04  | 2.442e+04  | 2.956e+04 | 
| TDIV1   | 6.270e+03  | 7.673e+03  | 6.269e+03  | 9.702e+03 | 
| NPOP2   | 6.819e+04  | 8.118e+04  | 6.934e+04  | 9.718e+04 | 
| NPOP1   | 1.666e+04  | 1.940e+04  | 1.736e+04  | 2.226e+04 | 
| MIG10R  | 1.175e-04  | 9.551e-05  | 7.034e-05  | 1.119e-04 | 
| MIG01R  | 9.502e-05  | 7.575e-05  | 5.929e-05  | 8.985e-05 | 
| MIG12R  | 2.868e-04  | 2.326e-04  | 1.710e-04  | 2.761e-04 | 
| MIG21R  | 8.777e-05  | 7.292e-05  | 6.057e-05  | 8.378e-05 | 
| MIG32R  | 6.628e-05  | 5.836e-05  | 5.131e-05  | 6.547e-05 | 
| MIG23R  | 8.825e-06  | 8.822e-06  | 8.757e-06  | 8.826e-06 | 
| TMRMG   | 5.924e+03  | 7.138e+03  | 6.051e+03  | 8.578e+03 | 
| recomb  | 3.204e-09  | 3.210e-09  | 3.071e-09  | 3.229e-09 | 
| TEGR2   | 1.266e+04  | 2.083e+04  | 1.170e+04  | 3.655e+04 | 
| TEGR1   | 3.500e+03  | 3.874e+03  | 8.631e+02  | 9.964e+03 | 
| GROW1   | -4.940e-05 | -6.776e-05 | -1.990e-05 | -9.709e-05 | 

The bootstrap corrections and confidence ranges presented in the
table are based on finite-sites bootstrap simulations.

---

# Argweaver

<details open>
<summary>argweaver</summary>
<br>
Prepare the site input file needed for argweaver. We selected scaff 78, 60 and 70 for chr18. I include 10 "neutral" scaffolds: scaff68 from chr15, scaff120, 174, 93 from chr21, scaff83, 125, 295 from chr19, and scaff103, 123, 129 from chr20. These chromosome have similar GC as chr 18 (47-49%)

Run `1.9.0_pre_argweaver.sh`

```
#!/bin/bash -l
#SBATCH -J preargweaver
#SBATCH --get-user-env
#SBATCH --mail-user=gwee@biologie.uni-muenchen.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=2-00:00:00
#SBATCH -o /dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2/slurms/slurm-%
j-%x.out

# sbatch 1.9.0_pre_argweaver.sh

module load vcftools/0.1.14-gcc8
module load bcftools/1.10.2-gcc8
conda activate biotools

ref="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2/genome_HC_allpaths41687_v2.5_chrW.fasta"
dat="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2"

echo $(date)
STARTTIME=$(date +%s)

cd ${dat}/argweaver

#create a clean geno file with invariant sites included. Missing sites due to repeat regions
#subset vcf to 4 pops (germany, spain, IRQ, hooded)
vcftools --keep poplist.txt --gzvcf ${dat}/05.1_recal/overlap/134inds_overlapped_filtered_norepeats_hwe_incinvariants_norepeats.vcf.gz --recode --out 4pop_hwe_incinvariants_norepeats
bcftools view -Oz -o 4pop_hwe_incinvariants_norepeats.vcf.gz 4pop_hwe_incinvariants_norepeats.recode.vcf
bcftools index 4pop_hwe_incinvariants_norepeats.vcf.gz
rm 4pop_hwe_incinvariants_norepeats.recode.vcf

python ${dat}/06_results/fst/genomics_general/VCF_processing/parseVCF.py -i 4pop_hwe_incinvariants_norepeats.vcf.gz --skipIndels --ploidyMismatchToMissing --ploidy 2 -o 4pop_hwe_incinvariants_norepeats.geno

sed -e 's////g' 4pop_hwe_incinvariants_norepeats.geno > 4pop_hwe_incinvariants_norepeats_edited.geno

for i in 120 174 93 73 125 103 123 129 295 68 60 70 78
do
grep -w "scaffold_${i}" 4pop_hwe_incinvariants_norepeats_edited.geno > 4pop_hwe_incinvariants_norepe
ats.geno.scaff${i}
awk '{print $2}' 4pop_hwe_incinvariants_norepeats.geno.scaff${i} > 4pop_hwe_incinvariants_norepeats.
geno.scaff${i}.sites
done

#germany(15; col3-17), spain(15; col18-32), IRQ(5; col33-37), hooded(15; col38-52)

#for cor1 and cnx3
for i in 120 174 93 73 125 103 123 129 295 68 60 70 78
do
mkdir scaff${i}_cor1cnx3
cut -d$'\t' -f18-22,38-42 4pop_hwe_incinvariants_norepeats.geno.scaff${i} | sed -e 's/\t//g' | sed -
e 's/\///g' > 4pop_hwe_incinvariants_norepeats.geno.scaff${i}.cor1cnx3.tmp
awk '
    {
        for ( i = 1; i <= NF; i++ )
            {
                if ( length ($i) > 20 )
                    {
                        gsub($i, "NNNNNNNNNNNNNNNNNNNN", $0)
                    }
            }
    }
1 ' 4pop_hwe_incinvariants_norepeats.geno.scaff${i}.cor1cnx3.tmp > 4pop_hwe_incinvariants_norepeats.
geno.scaff${i}.cor1cnx3.3.tmp
paste -d"\t" 4pop_hwe_incinvariants_norepeats.geno.scaff${i}.sites 4pop_hwe_incinvariants_norepeats.
geno.scaff${i}.cor1cnx3.3.tmp > 4pop_hwe_incinvariants_norepeats.geno.scaff${i}.cor1cnx3.2.tmp
cat ./header/header_cor1cnx3_scaff${i}.txt 4pop_hwe_incinvariants_norepeats.geno.scaff${i}.cor1cnx3.
2.tmp > ./scaff${i}_cor1cnx3/4pop_hwe_incinvariants_norepeats.geno.scaff${i}.cor1cnx3.sites
rm *cor1cnx3*.tmp
done

#for cnx3 and cnx6
for i in 120 174 93 73 125 103 123 129 295 68 60 70 78
do
mkdir scaff${i}_cnx3cnx6
cut -d$'\t' -f38-42,33-37 4pop_hwe_incinvariants_norepeats.geno.scaff${i} | sed -e 's/\t//g' | sed -
e 's/\///g' > 4pop_hwe_incinvariants_norepeats.geno.scaff${i}.cnx3cnx6.tmp
awk '
    {
        for ( i = 1; i <= NF; i++ )
            {
                if ( length ($i) > 20 )
                    {
                        gsub($i, "NNNNNNNNNNNNNNNNNNNN", $0)
                    }
            }
    }
1 ' 4pop_hwe_incinvariants_norepeats.geno.scaff${i}.cnx3cnx6.tmp > 4pop_hwe_incinvariants_norepeats.
geno.scaff${i}.cnx3cnx6.3.tmp
paste -d"\t" 4pop_hwe_incinvariants_norepeats.geno.scaff${i}.sites 4pop_hwe_incinvariants_norepeats.
geno.scaff${i}.cnx3cnx6.3.tmp > 4pop_hwe_incinvariants_norepeats.geno.scaff${i}.cnx3cnx6.2.tmp
cat ./header/header_cnx3cnx6_scaff${i}.txt 4pop_hwe_incinvariants_norepeats.geno.scaff${i}.cnx3cnx6.
2.tmp > ./scaff${i}_cnx3cnx6/4pop_hwe_incinvariants_norepeats.geno.scaff${i}.cnx3cnx6.sites
rm *cnx3cnx6*.tmp
done

#for cor1 and cor2
for i in 120 174 93 73 125 103 123 129 295 68 60 70 78
do
mkdir scaff${i}_cor1cor2
cut -d$'\t' -f3-7,18-22 4pop_hwe_incinvariants_norepeats.geno.scaff${i} | sed -e 's/\t//g' | sed -e 
's/\///g' > 4pop_hwe_incinvariants_norepeats.geno.scaff${i}.cor1cor2.tmp
awk '
    {
        for ( i = 1; i <= NF; i++ )
            {
                if ( length ($i) > 20 )
                    {
                        gsub($i, "NNNNNNNNNNNNNNNNNNNN", $0)
                    }
            }
    }
1 ' 4pop_hwe_incinvariants_norepeats.geno.scaff${i}.cor1cor2.tmp > 4pop_hwe_incinvariants_norepeats.
geno.scaff${i}.cor1cor2.3.tmp
paste -d"\t" 4pop_hwe_incinvariants_norepeats.geno.scaff${i}.sites 4pop_hwe_incinvariants_norepeats.
geno.scaff${i}.cor1cor2.3.tmp > 4pop_hwe_incinvariants_norepeats.geno.scaff${i}.cor1cor2.2.tmp
cat ./header/header_cor1cor2_scaff${i}.txt 4pop_hwe_incinvariants_norepeats.geno.scaff${i}.cor1cor2.
2.tmp > ./scaff${i}_cor1cor2/4pop_hwe_incinvariants_norepeats.geno.scaff${i}.cor1cor2.sites 
rm *cor1cor2*.tmp
done

#for cor1 and cnx6
for i in 120 174 93 73 125 103 123 129 295 68 60 70 78
do
mkdir scaff${i}_cor1cnx6
cut -d$'\t' -f18-22,33-37 4pop_hwe_incinvariants_norepeats.geno.scaff${i} | sed -e 's/\t//g' | sed -
e 's/\///g' > 4pop_hwe_incinvariants_norepeats.geno.scaff${i}.cor1cnx6.tmp
awk '
    {
        for ( i = 1; i <= NF; i++ )
            {
                if ( length ($i) > 20 )
                    {
                        gsub($i, "NNNNNNNNNNNNNNNNNNNN", $0)
                    }
            }
    }
1 ' 4pop_hwe_incinvariants_norepeats.geno.scaff${i}.cor1cnx6.tmp > 4pop_hwe_incinvariants_norepeats.
geno.scaff${i}.cor1cnx6.3.tmp
paste -d"\t" 4pop_hwe_incinvariants_norepeats.geno.scaff${i}.sites 4pop_hwe_incinvariants_norepeats.
geno.scaff${i}.cor1cnx6.3.tmp > 4pop_hwe_incinvariants_norepeats.geno.scaff${i}.cor1cnx6.2.tmp
cat ./header/header_cor1cnx6_scaff${i}.txt 4pop_hwe_incinvariants_norepeats.geno.scaff${i}.cor1cnx6.
2.tmp > ./scaff${i}_cor1cnx6/4pop_hwe_incinvariants_norepeats.geno.scaff${i}.cor1cnx6.sites
rm *cor1cnx6*.tmp
done

#for cor2 and cnx3
for i in 120 174 93 73 125 103 123 129 295 68 60 70 78
do
mkdir scaff${i}_cor2cnx3
cut -d$'\t' -f3-7,38-42 4pop_hwe_incinvariants_norepeats.geno.scaff${i} | sed -e 's/\t//g' | sed -e 
's/\///g' > 4pop_hwe_incinvariants_norepeats.geno.scaff${i}.cor2cnx3.tmp
awk '
    {
        for ( i = 1; i <= NF; i++ )
            {
                if ( length ($i) > 20 )
                    {
                        gsub($i, "NNNNNNNNNNNNNNNNNNNN", $0)
                    }
            }
    }
1 ' 4pop_hwe_incinvariants_norepeats.geno.scaff${i}.cor2cnx3.tmp > 4pop_hwe_incinvariants_norepeats.
geno.scaff${i}.cor2cnx3.3.tmp
paste -d"\t" 4pop_hwe_incinvariants_norepeats.geno.scaff${i}.sites 4pop_hwe_incinvariants_norepeats.
geno.scaff${i}.cor2cnx3.3.tmp > 4pop_hwe_incinvariants_norepeats.geno.scaff${i}.cor2cnx3.2.tmp
cat ./header/header_cor2cnx3_scaff${i}.txt 4pop_hwe_incinvariants_norepeats.geno.scaff${i}.cor2cnx3.
2.tmp > ./scaff${i}_cor2cnx3/4pop_hwe_incinvariants_norepeats.geno.scaff${i}.cor2cnx3.sites
rm *cor2cnx3*.tmp
done

for i in 120 174 93 73 125 103 123 129 295 68 60 70 78
do
for pair in cor2cnx3 cor2cnx6 cor1cnx6 cor1cor2 cnx3cnx6 cor1cnx3
do
cd scaff${i}_${pair}
grep -v "TTTTTTTTTTTTTTTTTTTT\|CCCCCCCCCCCCCCCCCCCC\|AAAAAAAAAAAAAAAAAAAA\|GGGGGGGGGGGGGGGGGGGG" 4po
p_hwe_incinvariants_norepeats.geno.scaff${i}.${pair}.sites > 4pop_hwe_incinvariants_norepeats.geno.s
caff${i}.${pair}.variantonly.sites
cd ../
done
done

for i in 120 174 93 73 125 103 123 129 295 68 60 70 78
do
grep scaffold_${i} ${dat}/ref2.5repeats_header.bed > repeatmasked.scaff${i}.bed
sed -i "s/scaffold_/scaff/g" repeatmasked.scaff${i}.bed
done

ENDTIME=$(date +%s)
echo $(date)
echo "It takes $(($ENDTIME - $STARTTIME)) seconds to complete this task"
echo "It takes $((($ENDTIME - $STARTTIME)/60)) mins to complete this task"
echo "It takes $((($ENDTIME - $STARTTIME)/3600)) hours to complete this task"
```

Run `1.9.1_argweaver_unphased_masked.sh` for the neutral scaffolds and chr18 scaffolds
 ```
 #!/bin/bash -l
#SBATCH --get-user-env
#SBATCH --mail-user=gwee@biologie.uni-muenchen.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_production
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=7-00:00:00
#SBATCH -o /dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2/slurms/slurm-%
j-%x.out

# 6 pop pair comparisons
# for i in in scaff*_cor1cor2; do base=${i%_cor1cor2*}; sbatch --job-name=arg${i} /dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2/1.9.1_argweaver_unphased_masked.sh $i ${base}; done
# for i in in scaff*_cor2cnx6; do base=${i%_cor2cnx6*}; sbatch --job-name=arg${i} /dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2/1.9.1_argweaver_unphased_masked.sh $i ${base}; done
# for i in in scaff*_cor2cnx3; do base=${i%_cor2cnx3*}; sbatch --job-name=arg${i} /dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2/1.9.1_argweaver_unphased_masked.sh $i ${base}; done
##
# for i in in scaff*_cnx3cnx6; do base=${i%_cnx3cnx6*}; sbatch --job-name=arg${i} /dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2/1.9.1_argweaver_unphased_masked.sh $i ${base}; done
# for i in in scaff*_cor1cnx6; do base=${i%_cor1cnx6*}; sbatch --job-name=arg${i} /dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2/1.9.1_argweaver_unphased_masked.sh $i ${base}; done
# for i in in scaff*_cor1cnx3; do base=${i%_cor1cnx3*}; sbatch --job-name=arg${i} /dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2/1.9.1_argweaver_unphased_masked.sh $i ${base}; done

conda activate biotools

ref="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2/genome_HC_allpaths41687_v2.5_chrW.fasta"
dat="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2"

echo $(date)
STARTTIME=$(date +%s)

cd ${dat}/argweaver

cd $1
echo $1

arg-sample -s 4pop_hwe_incinvariants_norepeats.geno.${1}.variantonly.sites -o 4pop_hwe_incinvariants
_norepeats.geno.${1} -N 10000 \
-r 1.6e-9 -m 3.18e-9 --ntimes 20 --maxtime 1e7 -n 2000 --unphased --delta 0.005 \
--randomize-phase 1.0 --sample-phase 10 --resample-window-iters 1 --resample-window 10000 --maskmap 
../repeatmasked.${2}.bed

mkdir burnin
for i in {0..1000..10}
do 
mv 4pop_hwe_incinvariants_norepeats.geno.${1}.${i}.smc.gz ./burnin
mv 4pop_hwe_incinvariants_norepeats.geno.${1}.${i}.sites.gz ./burnin
done

# cd /home/gwee/argweaver
# for i in scaff*
# do
# python2 ./bin/arg-extract-tmrca ${i}/4pop_hwe_incinvariants_norepeats.geno.${i}.%d.smc.gz > ${i}.tmrca.txt
# done

ENDTIME=$(date +%s)
echo $(date)
echo "It takes $(($ENDTIME - $STARTTIME)) seconds to complete this task"
echo "It takes $((($ENDTIME - $STARTTIME)/60)) mins to complete this task"
echo "It takes $((($ENDTIME - $STARTTIME)/3600)) hours to complete this task"
```

Run `1.9.1_argweaver_unphased_masked_region.sh` for the peak region in scaffold 60 and 78

```
#!/bin/bash -l
#SBATCH -J argelevated
#SBATCH --get-user-env
#SBATCH --mail-user=gwee@biologie.uni-muenchen.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=1-00:00:00
#SBATCH -o /dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2/slurms/slurm-%
j-%x.out

# sbatch 1.9.1_argweaver_unphased_masked_region.sh cnx3cnx6
# sbatch 1.9.1_argweaver_unphased_masked_region.sh cor2cnx3
# sbatch 1.9.1_argweaver_unphased_masked_region.sh cor2cnx6
# sbatch 1.9.1_argweaver_unphased_masked_region.sh cor1cnx6
# sbatch 1.9.1_argweaver_unphased_masked_region.sh cor1cnx3
# sbatch 1.9.1_argweaver_unphased_masked_region.sh cor1cor2

conda activate biotools

ref="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2/genome_HC_allpaths41687_v2.5_chrW.fasta"
dat="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2"

echo $(date)
STARTTIME=$(date +%s)

cd ${dat}/argweaver

mkdir scaff78_${1}_elevated
cd scaff78_${1}_elevated
arg-sample -s ../scaff78_${1}/4pop_hwe_incinvariants_norepeats.geno.scaff78_${1}.variantonly.sites -o 4pop_hwe_incinvariants_norepeats.geno.scaff78elevated_${1} -N 10000 \
-r 1.6e-9 -m 3.18e-9 --ntimes 20 --maxtime 1e7 -n 2000 --unphased --delta 0.005 --region 1363447-2400000 \
--randomize-phase 1.0 --sample-phase 10 --resample-window-iters 1 --resample-window 10000 --maskmap ../repeatmasked.scaff78.bed

mkdir burnin
for bs in {0..1000..10}
do 
mv 4pop_hwe_incinvariants_norepeats.geno.scaff78elevated_${1}.${bs}.smc.gz ./burnin
mv 4pop_hwe_incinvariants_norepeats.geno.scaff78elevated_${1}.${bs}.sites.gz ./burnin
done

cd ${dat}/argweaver
mkdir scaff60_${1}_elevated
cd scaff60_${1}_elevated
arg-sample -s ../scaff60_${1}/4pop_hwe_incinvariants_norepeats.geno.scaff60_${1}.variantonly.sites -o 4pop_hwe_incinvariants_norepeats.geno.scaff60elevated_${1} -N 10000 \
-r 1.6e-9 -m 3.18e-9 --ntimes 20 --maxtime 1e7 -n 2000 --unphased --delta 0.005 --region 50000-607
285 \
--randomize-phase 1.0 --sample-phase 10 --resample-window-iters 1 --resample-window 10000 --maskma
p ../repeatmasked.scaff60.bed

mkdir burnin
for bs in {0..1000..10}
do 
mv 4pop_hwe_incinvariants_norepeats.geno.scaff60elevated_${1}.${bs}.smc.gz ./burnin
mv 4pop_hwe_incinvariants_norepeats.geno.scaff60elevated_${1}.${bs}.sites.gz ./burnin
done

# cd /home/gwee/argweaver
# for i in scaff*
# do
# python2 ./bin/arg-extract-tmrca ${i}/4pop_hwe_incinvariants_norepeats.geno.${i}.%d.smc.gz > ${i}.tmrca.txt
# done

ENDTIME=$(date +%s)
echo $(date)
echo "It takes $(($ENDTIME - $STARTTIME)) seconds to complete this task"
echo "It takes $((($ENDTIME - $STARTTIME)/60)) mins to complete this task"
echo "It takes $((($ENDTIME - $STARTTIME)/3600)) hours to complete this task"
```

The `*.tmrca.txt` files are used for plotting. There are 6 columns in the file: chromosome, start, end, posterior mean TMRCA (generations), lower 2.5 percentile TMRCA, and upper 97.5 percentile TMRCA.

```
library(tidyverse)
library(zoo)
safe_colorblind_palette <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888", "#000000", "#E69F00", "#56B4E9")

# create a function to generate random walk, which kind of creates LD-like patterns 
random_walk <- function(n, start=0, step_scale=0.05, outlier_prob=0.99, outlier_scale=0.9) {
  walk <- numeric(n)
  walk[1] <- start
  
  for(i in 2:n) {
    if(runif(1) < outlier_prob) {
      # generate outlier
      step <- rnorm(1, mean=0, sd=outlier_scale)
    } else {
      # generate regular step
      step <- rnorm(1, mean=0, sd=step_scale)
    }
    
    walk[i] <- walk[i-1] + step
  }
  
  # rescale to [0, 1] interval
  walk <- (walk - min(walk)) / (max(walk) - min(walk))
  
  return(walk)
}

setwd("/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2/argweaver")

scaff <- c("103","120","123","125","129","174","295","68","70","73","93","60","60_elevated","78","78_elevated")
length(scaff)
pairs <- c("cnx3cnx6","cor2cnx6","cor2cnx3","cor1cor2","cor1cnx3","cor1cnx6")
length(pairs)

results <- data.frame(matrix(data=NA,nrow=length(scaff),ncol=length(pairs)))
colnames(results) <- pairs

for (i in 1:length(scaff)) {
  for (p in 1:length(pairs)) {
    data1 <- read.delim(paste("scaff",scaff[i],"_",pairs[p],".tmrca.txt", sep=""), header=FALSE, sep="\t")
    sample_df <- tibble(
      V1 = data1$V1,                            
      V2 = data1$V3,    
      V3 = data1$V4)
    
    sample_df = sample_df %>% arrange(V1,V2) %>%
      group_by(V1) %>%  
      mutate(RollingMean = rollmean(V3, k = 500, align = 'center', fill = NA))
    
    rolling_mean = sample_df %>% 
      filter(!is.na(RollingMean)) %>%
      group_by(V1) %>%
      filter(row_number() %% 500 == 1)
    
    results[i,p] <- mean(data1$V4, na.rm=TRUE) 
    
    pdf(file=paste("tmrca_500bprolls_scaff",scaff[i],"_",pairs[p],".pdf",sep=""))
    print(rolling_mean %>% 
            ggplot(aes(x=V2,y=RollingMean))+
            geom_line(col='blue',lty=1,lwd=2)+
            theme_bw() + 
            scale_y_continuous(
              breaks = scales::pretty_breaks(n = 4))+  # Create 4 breaks
            xlab("bp position") + ylab(paste("scaff",scaff[i],"_",pairs[p],"TMRCA (gen)",sep="_"))) 
    dev.off()
    
  }
}

######################## combined plot ########################
setwd("/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2/argweaver")
d <- c("scaff78_elevated") 
d1 <- c("cor1cor2","cor2cnx3", "cor2cnx6","cnx3cnx6")

  data1 <- read.delim(paste(d,"_",d1[1],".tmrca.txt", sep=""), header=FALSE, sep="\t")
  sample_df <- tibble(
    V1 = data1$V1,                            
    V2 = data1$V3,    
    V3 = data1$V4)
  sample_df = sample_df %>% arrange(V1,V2) %>% 
    group_by(V1) %>%  
    mutate(RollingMean = rollmean(V3, k = 500, align = 'center', fill = NA)) 
  rolling_mean = sample_df %>% 
    filter(!is.na(RollingMean)) %>%
    group_by(V1) %>%
    filter(row_number() %% 500 == 1)

    data2 <- read.delim(paste(d,"_",d1[2],".tmrca.txt", sep=""), header=FALSE, sep="\t")
  sample_df2 <- tibble(
    V1 = data2$V1,                            
    V2 = data2$V3,    
    V3 = data2$V4)
  sample_df2 = sample_df2 %>% arrange(V1,V2) %>% 
    group_by(V1) %>% 
    mutate(RollingMean = rollmean(V3, k = 500, align = 'center', fill=NA))  
  rolling_mean2 = sample_df2 %>% 
    filter(!is.na(RollingMean)) %>%
    group_by(V1) %>%
    filter(row_number() %% 500 == 1)
  
  data3 <- read.delim(paste(d,"_",d1[3],".tmrca.txt", sep=""), header=FALSE, sep="\t")
  sample_df3 <- tibble(
    V1 = data3$V1,                            
    V2 = data3$V3,    
    V3 = data3$V4)
  sample_df3 = sample_df3 %>% arrange(V1,V2) %>% 
    group_by(V1) %>% 
    mutate(RollingMean = rollmean(V3, k = 500, align = 'center', fill = NA)) 
  rolling_mean3 = sample_df3 %>% 
    filter(!is.na(RollingMean)) %>%
    group_by(V1) %>%
    filter(row_number() %% 500 == 1)

    data4 <- read.delim(paste(d,"_",d1[4],".tmrca.txt", sep=""), header=FALSE, sep="\t")
  sample_df4 <- tibble(
    V1 = data4$V1,                            
    V2 = data4$V3,    
    V3 = data4$V4)
  sample_df4 = sample_df4 %>% arrange(V1,V2) %>% 
    group_by(V1) %>% 
    mutate(RollingMean = rollmean(V3, k = 500, align = 'center', fill = NA)) 
  rolling_mean4 = sample_df4 %>% 
    filter(!is.na(RollingMean)) %>%
    group_by(V1) %>%
    filter(row_number() %% 500 == 1)

  p1 <- ggplot()+
  geom_line(data=rolling_mean, aes(x=rolling_mean$V2, y=rolling_mean$RollingMean, colour=d1[1]),lty=1,lwd=1)+
  geom_line(data=rolling_mean2, aes(x=rolling_mean2$V2, y=rolling_mean2$RollingMean, colour=d1[2]),lty=1,lwd=1) + 
  geom_line(data=rolling_mean3, aes(x=rolling_mean3$V2, y=rolling_mean3$RollingMean, colour=d1[3]),lty=2,lwd=1) +
  geom_line(data=rolling_mean4, aes(x=rolling_mean4$V2, y=rolling_mean4$RollingMean, colour=d1[4]),lty=2,lwd=1) + 
  scale_colour_manual(values=c("darkblue","black","darkgrey","#88CCEE"))+
  labs(colour = 'Population') +
  theme_bw() + guides(colour="none") +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 4)) + ylim(2.5e+04 ,  3e+05) +  
  scale_x_continuous(breaks = scales::pretty_breaks(n = 4)) +  
  xlab("scaff78_elevated_region (bp)") + ylab("TMRCA (gen)") 

d <- c("scaff60_elevated") 
d1 <- c("cor1cor2","cor2cnx3", "cor2cnx6","cnx3cnx6")
  data1 <- read.delim(paste(d,"_",d1[1],".tmrca.txt", sep=""), header=FALSE, sep="\t")
  sample_df <- tibble(
    V1 = data1$V1,                            
    V2 = data1$V3,    
    V3 = data1$V4)
  sample_df = sample_df %>% arrange(V1,V2) %>% #ensure its arranged by chr and pos 
    group_by(V1) %>% 
    mutate(RollingMean = rollmean(V3, k = 500, align = 'center', fill = NA)) 
  rolling_mean5 = sample_df %>% 
    filter(!is.na(RollingMean)) %>%
    group_by(V1) %>%
    filter(row_number() %% 500 == 1)
  
  data2 <- read.delim(paste(d,"_",d1[2],".tmrca.txt", sep=""), header=FALSE, sep="\t")
  sample_df2 <- tibble(
    V1 = data2$V1,                            
    V2 = data2$V3,    
    V3 = data2$V4)
  sample_df2 = sample_df2 %>% arrange(V1,V2) %>% 
    group_by(V1) %>% 
    mutate(RollingMean = rollmean(V3, k = 500, align = 'center', fill=NA))  
  rolling_mean6 = sample_df2 %>% 
    filter(!is.na(RollingMean)) %>%
    group_by(V1) %>%
    filter(row_number() %% 500 == 1)

    data3 <- read.delim(paste(d,"_",d1[3],".tmrca.txt", sep=""), header=FALSE, sep="\t")
  sample_df3 <- tibble(
    V1 = data3$V1,                            
    V2 = data3$V3,    
    V3 = data3$V4)
  sample_df3 = sample_df3 %>% arrange(V1,V2) %>% 
    group_by(V1) %>% 
    mutate(RollingMean = rollmean(V3, k = 500, align = 'center', fill = NA)) 
  rolling_mean7 = sample_df3 %>% 
    filter(!is.na(RollingMean)) %>%
    group_by(V1) %>%
    filter(row_number() %% 500 == 1)
  
  data4 <- read.delim(paste(d,"_",d1[4],".tmrca.txt", sep=""), header=FALSE, sep="\t")
  sample_df4 <- tibble(
    V1 = data4$V1,                            
    V2 = data4$V3,    
    V3 = data4$V4)
  sample_df4 = sample_df4 %>% arrange(V1,V2) %>% 
    group_by(V1) %>%  #ensure rolling mean is calculated by chr 
    mutate(RollingMean = rollmean(V3, k = 500, align = 'center', fill = NA)) 
  rolling_mean8 = sample_df4 %>% 
    filter(!is.na(RollingMean)) %>%
    group_by(V1) %>%
    filter(row_number() %% 500 == 1)
  
  p2 <- ggplot()+
    geom_line(data=rolling_mean5, aes(x=rolling_mean5$V2, y=rolling_mean5$RollingMean, colour=d1[1]),lty=1,lwd=1)+
    geom_line(data=rolling_mean6, aes(x=rolling_mean6$V2, y=rolling_mean6$RollingMean, colour=d1[2]),lty=1,lwd=1) + 
    geom_line(data=rolling_mean7, aes(x=rolling_mean7$V2, y=rolling_mean7$RollingMean, colour=d1[3]),lty=2,lwd=1) +
    geom_line(data=rolling_mean8, aes(x=rolling_mean8$V2, y=rolling_mean8$RollingMean, colour=d1[4]),lty=2,lwd=1) + 
    labs(colour = 'Population pair') +
    theme_bw() + 
    scale_colour_manual(values=c("darkblue","black","darkgrey","#88CCEE"))+
    scale_y_continuous(breaks = scales::pretty_breaks(n = 4)) + ylim(2.5e+04 ,  3e+05) + 
    scale_x_continuous(breaks = scales::pretty_breaks(n = 4)) +
    xlab("scaff60_elevated_region (bp)") + theme(axis.title.y=element_blank(),axis.text.y=element_blank())

  grid.arrange(b1, arrangeGrob(p1, p2,ncol=2), nrow=2)
```
</details>

