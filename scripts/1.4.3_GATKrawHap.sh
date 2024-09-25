#!/bin/bash -l
#SBATCH -J GATKrawHap
#SBATCH --get-user-env
#SBATCH --mail-user=gwee@biologie.uni-muenchen.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_production
#SBATCH --cpus-per-task=4
#SBATCH --time=7-00:00:00
#SBATCH -o /dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2/slurms/slurm-%j-%x.out

# cd 04_markdup
# for fname in *_markdup_cleaned.bam; do base=${fname%_markdup_cleaned.bam*}; sbatch ../1.4.3_GATKrawHap.sh ${base}; done
# for fname in  C.torquatus_X*_markdup_cleaned.bam; do base=${fname%_markdup_cleaned.bam*}; sbatch ../1.4.3_GATKrawHap.sh ${base}; done
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

#no indel realignment on GATK4
#not using Spark as it is still beta and not comparable, but has been tested and seem to run in parallel with `--spark-master local[*]`
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


# echo HaplotypeCaller
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

ENDTIME=$(date +%s)
echo $(date)
echo "It takes $(($ENDTIME - $STARTTIME)) seconds to complete this task"
echo "It takes $((($ENDTIME - $STARTTIME)/60)) mins to complete this task"
echo "It takes $((($ENDTIME - $STARTTIME)/3600)) hours to complete this task"
