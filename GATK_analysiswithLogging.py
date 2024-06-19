import os
import logging

class mtDNAAnalysis:
    def __init__(self, ref_file, known_sites_file, aligned_reads_dir, reads_dir, results_dir, data_dir, read1_file, read2_file):
        self.ref_file = ref_file
        self.known_sites_file = known_sites_file
        self.aligned_reads_dir = aligned_reads_dir
        self.reads_dir = reads_dir
        self.results_dir = results_dir
        self.data_dir = data_dir
        self.read1_file = read1_file
        self.read2_file = read2_file

        self.create_directories()
        self.setup_logger()

        self.fastqc_path = "/Users/bhupi/Downloads/FastQC/fastqc"
        self.bwa_path = "/Users/bhupi/Downloads/bwa-0.7.17/bwa"
        self.gatk_path = "java -Dsamjdk.use_async_io_read_samtools=false -Dsamjdk.use_async_io_write_samtools=true -Dsamjdk.use_async_io_write_tribble=false -Dsamjdk.compression_level=2 -jar /Users/bhupi/Downloads/gatk-4.5.0.0/gatk-package-4.5.0.0-local.jar"

    def create_directories(self):
        """
        Create the necessary directories based on the file paths provided.
        """
        for directory in [self.aligned_reads_dir, self.reads_dir, self.results_dir, self.data_dir]:
            os.makedirs(directory, exist_ok=True)

    def setup_logger(self):
        """
        Set up the logger for the class.
        """
        self.logger = logging.getLogger(__name__)
        self.logger.setLevel(logging.DEBUG)

        # Create a file handler
        file_handler = logging.FileHandler('mtdna_analysis.log')
        file_handler.setLevel(logging.DEBUG)

        # Create a console handler
        console_handler = logging.StreamHandler()
        console_handler.setLevel(logging.INFO)

        # Create a formatter and add it to the handlers
        formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
        file_handler.setFormatter(formatter)
        console_handler.setFormatter(formatter)

        # Add the handlers to the logger
        self.logger.addHandler(file_handler)
        self.logger.addHandler(console_handler)

    def run_fastqc(self):
        self.logger.info("STEP 1: QC - Run fastqc")
        os.system(f"{self.fastqc_path} {self.reads_dir}/{self.read1_file} -o {self.reads_dir}/")
        os.system(f"{self.fastqc_path} {self.reads_dir}/{self.read2_file} -o {self.reads_dir}/")

    def map_to_reference(self):
        self.logger.info("STEP 2: Map to reference using BWA-MEM")
        os.system(f"{self.bwa_path} index {self.ref_file}")
        os.system(f"{self.bwa_path} mem -t 4 -R \"@RG\\tID:SRR062634\\tPL:ILLUMINA\\tSM:SRR062634\" {self.ref_file} {self.reads_dir}/NIJ2020-01J-08J-49-1-A_S7_L001_R1_001.fastq {self.reads_dir}/NIJ2020-01J-08J-49-1-A_S7_L001_R2_001.fastq > {self.aligned_reads_dir}/NIJ2020.paired.sam")

    def mark_duplicates_and_sort(self):
        self.logger.info("STEP 3: Mark Duplicates and Sort - GATK4")
        os.system(f"{self.gatk_path} MarkDuplicatesSpark -I {self.aligned_reads_dir}/NIJ2020.paired.sam -O {self.aligned_reads_dir}/NIJ2020_sorted_dedup_reads.bam")

    def base_quality_recalibration(self):
        self.logger.info("STEP 4: Base quality recalibration")
        os.system(f"{self.gatk_path} BaseRecalibrator -I {self.aligned_reads_dir}/NIJ2020_sorted_dedup_reads.bam -R {self.ref_file} --known-sites {self.known_sites_file} -O {self.data_dir}/recal_data.table")
        os.system(f"{self.gatk_path} ApplyBQSR -I {self.aligned_reads_dir}/NIJ2020_sorted_dedup_reads.bam -R {self.ref_file} --bqsr-recal-file {self.data_dir}/recal_data.table -O {self.aligned_reads_dir}/NIJ2020_sorted_dedup_bqsr_reads.bam")

    def collect_metrics(self):
        self.logger.info("STEP 5: Collect Alignment & Insert Size Metrics")
        os.system(f"{self.gatk_path} CollectAlignmentSummaryMetrics R={self.ref_file} I={self.aligned_reads_dir}/NIJ2020_sorted_dedup_bqsr_reads.bam O={self.aligned_reads_dir}/alignment_metrics.txt")
        os.system(f"{self.gatk_path} CollectInsertSizeMetrics INPUT={self.aligned_reads_dir}/NIJ2020_sorted_dedup_bqsr_reads.bam OUTPUT={self.aligned_reads_dir}/insert_size_metrics.txt HISTOGRAM_FILE={self.aligned_reads_dir}/insert_size_histogram.pdf")

    def call_variants(self):
        self.logger.info("STEP 6: Call Variants - gatk haplotype caller")
        os.system(
            f"{self.gatk_path} HaplotypeCaller -R {self.ref_file} -I {self.aligned_reads_dir}/NIJ2020_sorted_dedup_bqsr_reads.bam -O {self.results_dir}/raw_variants.vcf")

        self.logger.info("STEP 6: separating out the SNP and INDELs")
        os.system(
            f"{self.gatk_path} SelectVariants -R {self.ref_file} -V {self.results_dir}/raw_variants.vcf --select-type SNP -O {self.results_dir}/raw_snps.vcf")

    def run_all_steps(self):
        self.run_fastqc()
        self.map_to_reference()
        self.mark_duplicates_and_sort()
        self.base_quality_recalibration()
        self.collect_metrics()
        self.call_variants()

if __name__=="__main__":
    # Set the paths to the required files and directories
    mtdna_analysis = mtDNAAnalysis(
        ref_file="/Users/bhupi/Desktop/mtDNA/NIJ_mtDNA/supporting_files/mtDNA/chrM.fa",
        known_sites_file="/Users/bhupi/Desktop/mtDNA/NIJ_mtDNA/supporting_files/mtDNA/MITOMAP_HMTDB_known_indels_chrM.vcf",
        aligned_reads_dir="/Users/bhupi/Desktop/mtDNA/NIJ_mtDNA/aligned_reads",
        reads_dir="/Users/bhupi/Desktop/mtDNA/NIJ_mtDNA/reads",
        results_dir="/Users/bhupi/Desktop/mtDNA/NIJ_mtDNA/results",
        data_dir="/Users/bhupi/Desktop/mtDNA/NIJ_mtDNA/data",
        read1_file = "NIJ2020-01J-08J-49-1-A_S7_L001_R1_001.fastq",
        read2_file = "NIJ2020-01J-08J-49-1-A_S7_L001_R1_001.fastq"
    )

    mtdna_analysis.fastqc_path = "/Users/bhupi/Downloads/FastQC/fastqc"
    mtdna_analysis.bwa_path = "/Users/bhupi/Downloads/bwa-0.7.17/bwa"
    mtdna_analysis.gatk_path = "java -Dsamjdk.use_async_io_read_samtools=false -Dsamjdk.use_async_io_write_samtools=true -Dsamjdk.use_async_io_write_tribble=false -Dsamjdk.compression_level=2 -jar /Users/bhupi/Downloads/gatk-4.5.0.0/gatk-package-4.5.0.0-local.jar"

    mtdna_analysis.run_all_steps()
