###!/bin/bash
###
### Script to call germline variants in a human WGS paired end reads 2 X 100bp
### Following GATK4 best practices workflow - https://gatk.broadinstitute.org/hc/en-us/articles/360035535932-Germline-short-variant-discovery-SNPs-Indels-
### This script is for demonstration purposes only
#
#
# ##download data
##echo "Skip this STEP if you already have the data"
##wget ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/phase3/data/HG00096/sequence_read/SRR062634_1.filt.fastq.gz
##wget ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/phase3/data/HG00096/sequence_read/SRR062634_2.filt.fastq.gz
#
#
#echo "Run Prep files..."
#
###Prep files (TO BE GENERATED ONLY ONCE)
#
####download reference files
#echo "Downloading the reference File. Need to download only once"
#
#wget -P ~/Desktop/demo/supporting_files/hg38/ https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
#gunzip ~/Desktop/demo/supporting_files/hg38/hg38.fa.gz
#
##index ref - .fai file before running haplotype caller
#samtools faidx ~/Desktop/demo/supporting_files/hg38/hg38.fa
# go to samtool folder and then run this commands
#samtools-1.20 % ./samtools faidx ~/Desktop/mtDNA/NIJ_mtDNA/supporting_files/mtDNA/chrM.fa
#
# ref dict - .dict file before running haplotype caller
#gatk CreateSequenceDictionary R=~/Desktop/demo/supporting_files/hg38/hg38.fa O=~/Desktop/demo/supporting_files/hg38/hg38.dict
# bhupi@BhupindersMBP2 gatk-4.5.0.0 % ./gatk CreateSequenceDictionary R=/Users/bhupi/Desktop/mtDNA/NIJ_mtDNA/supporting_files/mtDNA/chrM.fa  O=/Users/bhupi/Desktop/mtDNA/NIJ_mtDNA/supporting_files/mtDNA/chrM.dict
#
# download known sites files for BQSR from GATK resource bundle
#wget -P ~/Desktop/demo/supporting_files/hg38/ https://storage.googleapis›.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf
#wget -P ~/Desktop/demo/supporting_files/hg38/ https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.idx
#
#
#
# VARIANT CALLING STEPS
#
#
# directories
ref="/Users/bhupi/Desktop/mtDNA/NIJ_mtDNA/supporting_files/mtDNA/chrM.fa"
##known_sites="/Users/bhupi/Desktop/demo/supporting_files/hg38/Homo_sapiens_assembly38.dbsnp138.vcf"
known_sites="/Users/bhupi/Desktop/mtDNA/NIJ_mtDNA/supporting_files/mtDNA/MITOMAP_HMTDB_known_indels_chrM.vcf"
#known_sites="/Users/bhupi/Desktop/mtDNA/NIJ_mtDNA/supporting_files/mtDNA/MITOMAP_HMTDB_known_indels_chrRSRS.vcf"
aligned_reads="/Users/bhupi/Desktop/mtDNA/NIJ_mtDNA/aligned_reads"
reads="/Users/bhupi/Desktop/mtDNA/NIJ_mtDNA/reads"
results="/Users/bhupi/Desktop/mtDNA/NIJ_mtDNA/results"
data="/Users/bhupi/Desktop/mtDNA/NIJ_mtDNA/data"
#
#
#
#
# -------------------
# STEP 1: QC - Run fastqcy
# -------------------
#
echo "STEP 1: QC - Run fastqc"
#
/Users/bhupi/Downloads/FastQC/fastqc ${reads}/SRR062634_1.filt.fastq.gz -o ${reads}/
/Users/bhupi/Downloads/FastQC/fastqc ${reads}/SRR062634_2.filt.fastq.gz -o ${reads}/
#
# No trimming required, quality looks okay.
#
#
# --------------------------------------
# STEP 2: Map to reference using BWA-MEM
# --------------------------------------
#
echo "STEP 2: Map to reference using BWA-MEM"
##
## BWA index reference
/Users/bhupi/Downloads/bwa-0.7.17/bwa  index ${ref}
##
##
## BWA alignment
/Users/bhupi/Downloads/bwa-0.7.17/bwa mem -t 4 -R  "@RG\tID:SRR062634\tPL:ILLUMINA\tSM:SRR062634" ${ref} ${reads}/NIJ2020-01J-08J-49-1-A_S7_L001_R1_001.fastq ${reads}/NIJ2020-01J-08J-49-1-A_S7_L001_R2_001.fastq > ${aligned_reads}/NIJ2020.paired.sam
##
##
##
##
## -----------------------------------------
## STEP 3: Mark Duplicates and Sort - GATK4
## -----------------------------------------
##
echo "STEP 3: Mark Duplicates and Sort - GATK4"
##
java -Dsamjdk.use_async_io_read_samtools=false -Dsamjdk.use_async_io_write_samtools=true -Dsamjdk.use_async_io_write_tribble=false -Dsamjdk.compression_level=2 -jar /Users/bhupi/Downloads/gatk-4.5.0.0/gatk-package-4.5.0.0-local.jar  MarkDuplicatesSpark -I ${aligned_reads}/NIJ2020.paired.sam -O ${aligned_reads}/NIJ2020_sorted_dedup_reads.bam
echo "STEP 3: Mark Duplicates and Sort - GATK4"
java -Dsamjdk.use_async_io_read_samtools=false -Dsamjdk.use_async_io_write_samtools=true -Dsamjdk.use_async_io_write_tribble=false -Dsamjdk.compression_level=2 -jar /Users/bhupi/Downloads/gatk-4.5.0.0/gatk MarkDuplicatesSpark -I ${aligned_reads}/NIJ2020.paired.sam -O ${aligned_reads}/NIJ2020_sorted_dedup_reads.bam
echo "STEP 3: Mark Duplicates and Sort - GATK4"
java -Dsamjdk.use_async_io_read_samtools=false -Dsamjdk.use_async_io_write_samtools=true -Dsamjdk.use_async_io_write_tribble=false -Dsamjdk.compression_level=2 -jar /Users/bhupi/Downloads/gatk-4.5.0.0/gatk-package-4.5.0.0-local.jar MarkDuplicatesSpark -I ${aligned_reads}/NIJ2020.paired.sam -O ${aligned_reads}/NIJ2020_sorted_dedup_reads.bam
##
##
#
#
# ----------------------------------
# STEP 4: Base quality recalibration
# ----------------------------------
#
##
echo "STEP 4: Base quality recalibration"
### at this step need to run this command ./gatk IndexFeatureFile -I /Users/bhupi/Desktop/mtDNA/NIJ_mtDNA/supporting_files/mtDNA/MITOMAP_HMTDB_known_indels_chrM.vcf
#
java -Dsamjdk.use_async_io_read_samtools=false -Dsamjdk.use_async_io_write_samtools=true -Dsamjdk.use_async_io_write_tribble=false -Dsamjdk.compression_level=2 -jar /Users/bhupi/Downloads/gatk-4.5.0.0/gatk-package-4.5.0.0-local.jar BaseRecalibrator -I ${aligned_reads}/NIJ2020_sorted_dedup_reads.bam -R ${ref} --known-sites ${known_sites} -O ${data}/recal_data.table
#
echo "STEP 4.2: Apply the model to adjust the base quality scores"
#
#
java -Dsamjdk.use_async_io_read_samtools=false -Dsamjdk.use_async_io_write_samtools=true -Dsamjdk.use_async_io_write_tribble=false -Dsamjdk.compression_level=2 -jar /Users/bhupi/Downloads/gatk-4.5.0.0/gatk-package-4.5.0.0-local.jar ApplyBQSR -I ${aligned_reads}/NIJ2020_sorted_dedup_reads.bam -R ${ref} --bqsr-recal-file ${data}/recal_data.table -O ${aligned_reads}/NIJ2020_sorted_dedup_bqsr_reads.bam

##
# -----------------------------------------------
# STEP 5: Collect Alignment & Insert Size Metrics
# -----------------------------------------------
#

echo "STEP 5: Collect Alignment & Insert Size Metrics"

java -Dsamjdk.use_async_io_read_samtools=false -Dsamjdk.use_async_io_write_samtools=true -Dsamjdk.use_async_io_write_tribble=false -Dsamjdk.compression_level=2 -jar /Users/bhupi/Downloads/gatk-4.5.0.0/gatk-package-4.5.0.0-local.jar  CollectAlignmentSummaryMetrics R=${ref} I=${aligned_reads}/NIJ2020_sorted_dedup_bqsr_reads.bam O=${aligned_reads}/alignment_metrics.txt

java -Dsamjdk.use_async_io_read_samtools=false -Dsamjdk.use_async_io_write_samtools=true -Dsamjdk.use_async_io_write_tribble=false -Dsamjdk.compression_level=2 -jar /Users/bhupi/Downloads/gatk-4.5.0.0/gatk-package-4.5.0.0-local.jar  CollectInsertSizeMetrics INPUT=${aligned_reads}/NIJ2020_sorted_dedup_bqsr_reads.bam OUTPUT=${aligned_reads}/insert_size_metrics.txt HISTOGRAM_FILE=${aligned_reads}/insert_size_histogram.pdf
java -Dsamjdk.use_async_io_read_samtools=false -Dsamjdk.use_async_io_write_samtools=true -Dsamjdk.use_async_io_write_tribble=false -Dsamjdk.compression_level=2 -jar /Users/bhupi/Downloads/gatk-4.5.0.0/gatk-package-4.5.0.0-local.jar CollectAlignmentSummaryMetrics R=${ref} I=${aligned_reads}/NIJ2020_sorted_dedup_bqsr_reads.bam O=${aligned_reads}/alignment_metrics.txt

echo "STEP 5_2: Collect Alignment & Insert Size Metrics"
echo "This step require R studio "
java -Dsamjdk.use_async_io_read_samtools=false -Dsamjdk.use_async_io_write_samtools=true -Dsamjdk.use_async_io_write_tribble=false -Dsamjdk.compression_level=2 -jar /Users/bhupi/Downloads/gatk-4.5.0.0/gatk-package-4.5.0.0-local.jar CollectInsertSizeMetrics INPUT=${aligned_reads}/NIJ2020_sorted_dedup_bqsr_reads.bam OUTPUT=${aligned_reads}/insert_size_metrics.txt HISTOGRAM_FILE=${aligned_reads}/insert_size_histogram.pdf

##
#
#
# ----------------------------------------------
# STEP 6: Call Variants - gatk haplotype caller
# ----------------------------------------------
#
echo "STEP 6: Call Variants - gatk haplotype caller"

java -Dsamjdk.use_async_io_read_samtools=false -Dsamjdk.use_async_io_write_samtools=true -Dsamjdk.use_async_io_write_tribble=false -Dsamjdk.compression_level=2 -jar /Users/bhupi/Downloads/gatk-4.5.0.0/gatk-package-4.5.0.0-local.jar HaplotypeCaller -R ${ref} -I ${aligned_reads}/NIJ2020_sorted_dedup_bqsr_reads.bam -O ${results}/raw_variants.vcf



echo "STEP 6: seperating out the SNP and INDELs"
java -Dsamjdk.use_async_io_read_samtools=false -Dsamjdk.use_async_io_write_samtools=true -Dsamjdk.use_async_io_write_tribble=false -Dsamjdk.compression_level=2 -jar /Users/bhupi/Downloads/gatk-4.5.0.0/gatk-package-4.5.0.0-local.jar SelectVariants -R ${ref} -V ${results}/raw_variants.vcf --select-type SNP -O ${results}/raw_snps.vcf
java -Dsamjdk.use_async_io_read_samtools=false -Dsamjdk.use_async_io_write_samtools=true -Dsamjdk.use_async_io_write_tribble=false -Dsamjdk.compression_level=2 -jar /Users/bhupi/Downloads/gatk-4.5.0.0/gatk-package-4.5.0.0-local.jar SelectVariants -R ${ref} -V ${results}/raw_variants.vcf --select-type INDEL -O ${results}/raw_indels.vcf

#
#
#
#
#
