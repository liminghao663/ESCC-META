
#############
# please replace "localtion_to_xxx" as the real paths in your computer
# please change the memory(such as "-Xmx100G") and thread parameters according to your computer
#####################

fastqdump=/localtion_to_sratoolkit/sratoolkit.2.11.3-ubuntu64/bin/fastq-dump
fastp=/localtion_to_fastp/fastp
GATK=/localtion_to_GATK4/gatk-4.1.6.0/gatk
GENOME=/localtion_to_genome/GRCh38.p13.genome.fa
INDEX=${GENOME}
PICARD=/localtion_to_picard/picard.jar
SNP1=/location_to_snp_dir/dbsnp_146.hg38.vcf.gz
SNP2=/location_to_snp_dir/1000G_phase1.snps.high_confidence.hg38.vcf.gz
SNP3=/location_to_snp_dir/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
samtools=/location_to_samtools/samtools
bwa=/location_to_samtools/bwa


#SRA file to fastq file
# please cd to the dir of downloaded SRA files

for in_sra in SRR*
do
echo $in_sra
$fastqdump --split-3 --gzip $in_sra
done


#QC of fastq

for name0 in *_1.fastq.gz
do
name=${name0/_1.fastq.gz/}
in1=${name}_1.fastq.gz
in2=${name}_2.fastq.gz
out1=${name}_R1.fastq.gz
out2=${name}_R2.fastq.gz

$fastp -i $in1 -o $out1 -I $in2 -O $out2 -h ${name}.html --json ${name}.json --thread 16
done

############################
#Mapping to genome

mkdir res
out_dir=res

for fq_gz in *_R1.fastq.gz
do
sample=${fq_gz/_R1.fastq.gz/}
in_fq1=${sample}_R1.fastq.gz
in_fq2=${sample}_R2.fastq.gz

$bwa mem \
-t 60 \
-M \
-v 2 \
-R "@RG\tID:$sample\tSM:$sample\tLB:WGS\tPL:WGS" \
$INDEX \
$in_fq1 $in_fq2 \
>$out_dir/$sample.sam

java -Xmx100g -jar $PICARD SortSam \
SORT_ORDER coordinate \
INPUT $out_dir/$sample.sam \
OUTPUT $out_dir/$sample.bam

$samtools flagstat --threads 60 $out_dir/$sample.bam > $out_dir/$sample.bam.alignment.flagstat
$samtools stats --threads 60 $out_dir/$sample.bam > $out_dir/$sample.bam.alignment.stat

$GATK MarkDuplicates \
--INPUT $out_dir/$sample.bam \
--java-options "-Xmx100G" \
--OUTPUT $out_dir/$sample.deduped.bam \
--METRICS_FILE $out_dir/$sample.deduped.matirx.txt \
--REMOVE_DUPLICATES TRUE

$GATK BuildBamIndex \
--INPUT $out_dir/$sample.deduped.bam

$GATK BaseRecalibrator \
--input $out_dir/$sample.deduped.bam \
--java-options "-Xmx100G" \
--output $out_dir/$sample.recal_data.table \
--known-sites $SNP1 \
--known-sites $SNP2 \
--known-sites $SNP3 \
--reference $GENOME

$GATK ApplyBQSR \
--input $out_dir/$sample.deduped.bam \
--java-options "-Xmx100G" \
-bqsr $out_dir/$sample.recal_data.table \
--output $out_dir/$sample.recal.bam

rm $out_dir/$sample.deduped.bam
rm $out_dir/$sample.bam
rm $out_dir/$sample.sam

done

#####################
## call snv

cd $out_dir

function call_snv()
{

tumor_bams=${1}
tumor_bams_line=$(echo $tumor_bams |awk -F "," '{for(i=1;i<=NF;i++){out="-I " $i ".recal.bam ";printf out}}')

normal_bam=${2}.recal.bam
sample_id=${3}
normal_id=${2}

$GATK Mutect2 --java-options "-Xmx60G" -R $GENOME --native-pair-hmm-threads 25 $tumor_bams_line -I $normal_bam -normal $normal_id -O ${sample_id}.vcf

$GATK FilterMutectCalls \
-R $GENOME \
-V ${sample_id}.vcf \
-O ${sample_id}.filtered.vcf

$GATK SelectVariants -R $GENOME \
--exclude-filtered true \
-V ${sample_id}.filtered.vcf \
-O ${sample_id}.PASS.vcf 

}

### this function is applicable to both paired tumor-normal sample or multi- tumor samples with one normal control
# for paired tumor-normal sample: call_snv tumor normal id
# such as 
# callsnv SRR1056922 SRR1056907 BGI-ESCC-038T
# for multi- tumor samples with one normal control: call_snv tumor1,tumor2,tumor3,…… normal id
# such as 
# callsnv SRR7343743,SRR7343744,SRR7343745,SRR7343746 SRR7343751 ESCC8
