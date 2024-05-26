#!/bin/bash 

if false
then


samtools faidx /home/bharath/gatk/hg38/hg38.fa

gatk CreateSequenceDictionary R=/home/bharath/gatk/hg38/hg38.fa O=/home/bharath/gatk/hg38/hg38.dict

wget -P /home/bharath/gatk/hg38/ https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf

wget -P /home/bharath/gatk/hg38/ https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.idx


fi 



# directories
ref="/home/bharath/gatk/hg38/hg38.fa"
known_sites="/home/bharath/gatk/hg38/Homo_sapiens_assembly38.dbsnp138.vcf"
aligned_reads="/home/bharath/gatk/aligned_reads"
reads="/home/bharath/gatk/reads"
results="/home/bharath/gatk/results"
data="/home/bharath/gatk/data"


echo "STEP 1: QC - Run fastqc"

fastqc ${reads}/SRR062634_1.filt.fastq.gz -o ${reads}/
fastqc ${reads}/SRR062634_2.filt.fastq.gz -o ${reads}/

echo "STEP 2: Map to reference using BWA-MEM"

BWA index reference 
bwa index ${ref}


# BWA alignment
bwa mem -t 4 -R "@RG\tID:SRR062634\tPL:ILLUMINA\tSM:SRR062634" ${ref} ${reads}/SRR062634_1.filt.fastq.gz ${reads}/SRR062634_2.filt.fastq.gz > ${aligned_reads}/SRR062634.paired.sam

echo "STEP 3: Mark Duplicates and sort - gatk"

gatk MarkDuplicatesSpark -I ${aligned_reads}/SRR062634.paired.sam -O ${aligned_reads}/SRR062634_sorted_dedup_reads.bam

echo "Step 4: Base quality recalibration"

gatk BaseRecalibrator -I ${aligned_reads}/SRR062634_sorted_dedup_reads.bam -R ${ref} --known-sites ${known_sites} -O ${data}/recal_data.table

gatk ApplyBQSR -I ${aligned_reads}/SRR062634_sorted_dedup_reads.bam -R ${ref} --bqsr-recal-file {$data}/recal_data.table -O ${aligned_reads}/SRR062634_sorted_dedup_bqsr_reads.bam 

echo "STEP 5: Collect Alignment & Insert Size Metrics"

gatk CollectAlignmentSummaryMetrics R=${ref} I=${aligned_reads}/SRR062634_sorted_dedup_bqsr_reads.bam O=${aligned_reads}/alignment_metrics.txt

echo "STEP 6: Call Variants - gatk haplotypecaller"
gatk HaplotypeCaller -R ${ref} -I ${aligned_reads}/SRR062634_sorted_dedup_bqsr_reads.bam -O ${results}/raw_variants.vcf
