import os
import logging
import sys
import tkinter as tk
from tkinter import filedialog
import sys

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
        #file_handler.setLevel(logging.DEBUG)
        file_handler.setLevel(logging.DEBUG)

        # Create a console handler
        console_handler = logging.StreamHandler()
        console_handler.setLevel(logging.DEBUG)

        # Create a formatter and add it to the handlers
        formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
        file_handler.setFormatter(formatter)
        console_handler.setFormatter(formatter)

        # Add the handlers to the logger
        self.logger.addHandler(file_handler)
        self.logger.addHandler(console_handler)



    def run_fastqc(self):
        self.logger.info("STEP 1 of 7: QC - Run fastqc")
       
        os.system(f"{self.fastqc_path} {self.reads_dir}/{self.read1_file} -o {self.reads_dir}/")
        os.system(f"{self.fastqc_path} {self.reads_dir}/{self.read2_file} -o {self.reads_dir}/")

    def map_to_reference(self):
        self.logger.info("STEP 2 of 7: Map to reference using BWA-MEM")
       

        os.system(f"{self.bwa_path} index {self.ref_file}")
        os.system(f"{self.bwa_path} mem -t 4 -R \"@RG\\tID:SRR062634\\tPL:ILLUMINA\\tSM:SRR062634\" {self.ref_file} {self.reads_dir}/NIJ2020-01J-08J-49-1-A_S7_L001_R1_001.fastq {self.reads_dir}/NIJ2020-01J-08J-49-1-A_S7_L001_R2_001.fastq > {self.aligned_reads_dir}/NIJ2020.paired.sam")

    def mark_duplicates_and_sort(self):
        self.logger.info("STEP 3 of 7: Mark Duplicates and Sort - GATK4")
       
        os.system(f"{self.gatk_path} MarkDuplicatesSpark -I {self.aligned_reads_dir}/NIJ2020.paired.sam -O {self.aligned_reads_dir}/NIJ2020_sorted_dedup_reads.bam")

    def base_quality_recalibration(self):
        self.logger.info("STEP 4 of 7: Base quality recalibration")
       
        os.system(f"{self.gatk_path} BaseRecalibrator -I {self.aligned_reads_dir}/NIJ2020_sorted_dedup_reads.bam -R {self.ref_file} --known-sites {self.known_sites_file} -O {self.data_dir}/recal_data.table")
        os.system(f"{self.gatk_path} ApplyBQSR -I {self.aligned_reads_dir}/NIJ2020_sorted_dedup_reads.bam -R {self.ref_file} --bqsr-recal-file {self.data_dir}/recal_data.table -O {self.aligned_reads_dir}/NIJ2020_sorted_dedup_bqsr_reads.bam")

    def collect_metrics(self):
        self.logger.info("STEP 5 of 7: Collect Alignment & Insert Size Metrics")
       
        os.system(f"{self.gatk_path} CollectAlignmentSummaryMetrics R={self.ref_file} I={self.aligned_reads_dir}/NIJ2020_sorted_dedup_bqsr_reads.bam O={self.aligned_reads_dir}/alignment_metrics.txt")
        os.system(f"{self.gatk_path} CollectInsertSizeMetrics INPUT={self.aligned_reads_dir}/NIJ2020_sorted_dedup_bqsr_reads.bam OUTPUT={self.aligned_reads_dir}/insert_size_metrics.txt HISTOGRAM_FILE={self.aligned_reads_dir}/insert_size_histogram.pdf")

    def call_variants(self):
        self.logger.info("STEP 6 of 7: Call Variants - gatk haplotype caller")
        os.system(f"{self.gatk_path} HaplotypeCaller -R {self.ref_file} -I {self.aligned_reads_dir}/NIJ2020_sorted_dedup_bqsr_reads.bam -O {self.results_dir}/raw_variants.vcf")

        self.logger.info("STEP 7 of 7: separating out the SNP and INDELs")
        os.system(f"{self.gatk_path} SelectVariants -R {self.ref_file} -V {self.results_dir}/raw_variants.vcf --select-type SNP -O {self.results_dir}/raw_snps.vcf")
        os.system(f"{self.gatk_path} SelectVariants -R {self.ref_file} -V {self.results_dir}/raw_variants.vcf --select-type INDEL -O {self.results_dir}/raw_indels.vcf")

    def run_all_steps(self):
        self.run_fastqc()
        self.map_to_reference()
        self.mark_duplicates_and_sort()
        self.base_quality_recalibration()
        self.collect_metrics()
        self.call_variants()


def choose_file():
    return filedialog.askopenfilename()


def choose_directory():
    return filedialog.askdirectory()


def main():
    root = tk.Tk()
    root.withdraw()  # Hide the root window


    # Ask user to select files and directories
    ref_file = choose_file()
    known_sites_file = choose_file()
    read1_file = choose_file()
    read2_file = choose_file()

    base_path = filedialog.askdirectory(title="Select Base Directory")

    aligned_reads_dir = os.path.join(base_path, "aligned_reads")
    reads_dir = os.path.join(base_path, "reads")
    results_dir = os.path.join(base_path, "results")
    data_dir = os.path.join(base_path, "data")

    mtdna_analysis = mtDNAAnalysis(
        ref_file=ref_file,
        known_sites_file=known_sites_file,
        aligned_reads_dir=aligned_reads_dir,
        reads_dir=reads_dir,
        results_dir=results_dir,
        data_dir=data_dir,
        read1_file=os.path.basename(read1_file),
        read2_file=os.path.basename(read2_file)
    )

    mtdna_analysis.run_all_steps()


if __name__ == "__main__":
    # with open("debug.log", "w") as f:
    #     # Redirect stdout to the file
    #     sys.stdout = f

        # Your code here
    main()

        # Restore stdout to its original value
        #sys.stdout = sys.__stdout__

