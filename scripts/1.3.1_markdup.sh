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

#save space
rm $1_markdup.txt
rm $1_markdup.bam
cd ${dat}/03_readgrp
rm *.bam


ENDTIME=$(date +%s)
echo $(date)
echo "It takes $(($ENDTIME - $STARTTIME)) sec to complete this task"
HOUR=3600
echo "It takes $((($ENDTIME - $STARTTIME) / $HOUR)) hr to complete this task"
