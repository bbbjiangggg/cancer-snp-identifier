import os
import time
import pyfiglet
from termcolor import colored
from snp_analysis_pipeline_v2.file_handling import read_accession_numbers, detect_accession_list_file
from snp_analysis_pipeline_v2.data_processing import detect_read_type, prefetch_and_convert, trim_reads
from snp_analysis_pipeline_v2.command_execution import run_command
from snp_analysis_pipeline_v2.path_management import ensure_directory, print_chromosome_paths
from snp_analysis_pipeline_v2.logging_module import log_message

# Function to delete intermediate files
def delete_intermediate_files(accession_number, chromosome):
    """Delete intermediate files like SAM, BAM, and raw BCF files."""
    intermediate_files = [
        f"{accession_number}/{accession_number}_mapped_{chromosome}.sam",
        f"{accession_number}/{accession_number}_mapped_{chromosome}.bam",
        f"{accession_number}/{accession_number}_mapped_{chromosome}.raw.bcf"
    ]
    
    for file_path in intermediate_files:
        if os.path.isfile(file_path):
            os.remove(file_path)
            log_message(f"Deleted {file_path}", level="success")
        else:
            log_message(f"File {file_path} not found. Skipping deletion.", level="warning")

def print_banner():
    banner_text = pyfiglet.figlet_format("CANCER IMMUNOLOGY", font="slant")
    print(colored(banner_text, "white"))
    sub_text = "SNP Analysis Pipeline v2.0"
    print(colored(f"{'='*len(sub_text)}", "yellow"))
    print(colored(sub_text.center(len(sub_text)), "green", attrs=["bold"]))
    print(colored(f"{'='*len(sub_text)}", "yellow"))

def main():
    # Print banner
    print_banner()

    bwa_base_path = "/usr/local/bin/bwa/"
    bowtie_base_path = "/usr/local/bin/bowtie/"
    fastp_path = "/usr/bin/fastp"

    # Step 1: Detect and select the accession list file (.txt)
    accession_list_file = detect_accession_list_file()  # This function already handles the print messages

    print()  # Line space after step 1

    # Automatically use the recommended number of threads
    available_threads = os.cpu_count()
    recommended_threads = max(1, available_threads // 2)
    log_message(f"Using {recommended_threads} threads for the analysis.", level="info")
    threads = str(recommended_threads)
    
    print()  # Line space after thread information

    # Step 2: Enter chromosomes to analyze
    print("2. Enter chromosomes to analyze (comma-separated) or 'all' for all chromosomes:", end=" ")
    chromosomes_input = input().strip()
    chromosomes_list = [str(i) for i in range(1, 23)] + ['X', 'Y'] if chromosomes_input.lower() == 'all' else chromosomes_input.split(',')

    print()  # Line space after chromosome input

    accession_numbers = read_accession_numbers(accession_list_file)

    # Check which accession numbers require analysis
    incomplete_accessions = []
    for accession_number in accession_numbers:
        analysis_needed = False
        for chromosome in chromosomes_list:
            final_vcf_file = f"{accession_number}/{accession_number}_mapped_{chromosome}.var.-final.vcf"
            if not os.path.isfile(final_vcf_file) or os.path.getsize(final_vcf_file) == 0:
                analysis_needed = True
                break
        if analysis_needed:
            incomplete_accessions.append(accession_number)

    total_incomplete = len(incomplete_accessions)
    log_message(f"Total incomplete accession numbers found: {total_incomplete}", level="info")
    
    print()  # Line space after showing total incomplete accessions

    # Ask the user how many accession numbers they wish to analyze
    if total_incomplete > 0:
        print(f"3. How many accession numbers would you like to analyze? (1-{total_incomplete}):", end=" ")
        num_to_analyze = int(input().strip())
        accession_numbers_to_analyze = incomplete_accessions[:num_to_analyze]
    else:
        log_message("All accession numbers are complete. Nothing to analyze.", level="success")
        return

    print()  # Line space after choosing how many accessions to analyze

    for accession_number in accession_numbers_to_analyze:
        ensure_directory(accession_number)

        # Detect if trimmed files already exist
        trimmed_file_single = f"{accession_number}/{accession_number}_trimmed.fq"
        trimmed_file_paired_1 = f"{accession_number}/{accession_number}_1_trimmed.fq"
        trimmed_file_paired_2 = f"{accession_number}/{accession_number}_2_trimmed.fq"

        if os.path.isfile(trimmed_file_single) or (os.path.isfile(trimmed_file_paired_1) and os.path.isfile(trimmed_file_paired_2)):
            log_message(f"Trimmed files already exist for {accession_number}. Skipping download and trimming.", level="warning")
            read_type = detect_read_type(accession_number)
        else:
            prefetch_and_convert(accession_number, threads)
            read_type = detect_read_type(accession_number)
            trim_reads(accession_number, read_type, fastp_path)

        # Line space after prefetching or skipping trimming
        print()

        # Perform the bowtie2 and bwa commands for the analysis, even if trimmed files exist
        for chromosome in chromosomes_list:
            log_message(f"Processing chromosome {chromosome} for {accession_number}...", level="info")

            # Start the timer
            start_time = time.time()

            input_file_for_mapping = trimmed_file_single if read_type == '1' else trimmed_file_paired_1

            log_message(f"\nMapping {accession_number} reads using Bowtie2 for chromosome {chromosome}...", level="info")
            bowtie_command = f"bowtie2 --very-fast-local -x {bowtie_base_path}hg38_index -U {input_file_for_mapping} -S {accession_number}/{accession_number}_mapped_{chromosome}.sam"
            run_command(bowtie_command)

            # Line space after mapping step
            print()

            # Convert SAM to BAM
            log_message(f"Converting SAM to BAM for chromosome {chromosome}...", level="info")
            run_command(f"samtools view -S -b {accession_number}/{accession_number}_mapped_{chromosome}.sam > {accession_number}/{accession_number}_mapped_{chromosome}.bam")

            # Line space after SAM to BAM conversion
            print()

            # Sort BAM file
            log_message(f"Sorting BAM for chromosome {chromosome}...", level="blue")
            run_command(f"samtools sort {accession_number}/{accession_number}_mapped_{chromosome}.bam > {accession_number}/{accession_number}_mapped_{chromosome}.sorted.bam")

            # Line space after BAM sorting
            print()

            # Generate mpileup and VCF
            log_message(f"Generating VCF for chromosome {chromosome}...", level="blue")
            bwa_chrom_path = f"/usr/local/bin/bwa/hg38.fa"
            run_command(f"bcftools mpileup -f {bwa_chrom_path} {accession_number}/{accession_number}_mapped_{chromosome}.sorted.bam | bcftools call -mv -Ob -o {accession_number}/{accession_number}_mapped_{chromosome}.raw.bcf")

            # Line space after generating VCF
            print()

            # Finalize VCF file
            log_message(f"Finalizing VCF for chromosome {chromosome}...", level="magenta")
            final_vcf_file = f"{accession_number}/{accession_number}_mapped_{chromosome}.var.-final.vcf"
            run_command(f"bcftools view {accession_number}/{accession_number}_mapped_{chromosome}.raw.bcf | vcfutils.pl varFilter - > {final_vcf_file}")

            # Line space after finalizing VCF
            print()

            # Delete intermediate files
            delete_intermediate_files(accession_number, chromosome)

            # End the timer and calculate elapsed time
            end_time = time.time()
            elapsed_time = end_time - start_time
            log_message(f"Time elapsed for chromosome {chromosome}: {elapsed_time:.2f} seconds", level="info")

            # Line space after each chromosome processing
            print()

if __name__ == "__main__":
    main()
