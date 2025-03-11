import os
import time
import pyfiglet
import subprocess
from termcolor import colored
from tqdm import tqdm
from snp_analysis_pipeline_v3.file_handling import read_accession_numbers, detect_accession_list_file
from snp_analysis_pipeline_v3.data_processing import detect_read_type, prefetch_and_convert, trim_reads
from snp_analysis_pipeline_v3.command_execution import run_command
from snp_analysis_pipeline_v3.path_management import ensure_directory, print_chromosome_paths
from snp_analysis_pipeline_v3.logging_module import log_message

# Function to delete intermediate files
def delete_intermediate_files(accession_number, chromosome):
    files_to_delete = [
        f"{accession_number}/{accession_number}_mapped_{chromosome}.sam",
        f"{accession_number}/{accession_number}_mapped_{chromosome}.bam",
        f"{accession_number}/{accession_number}_mapped_{chromosome}.raw.bcf",
        f"{accession_number}/{accession_number}_1.fastq",
        f"{accession_number}/{accession_number}_2.fastq",
        f"{accession_number}/{accession_number}.fastq",
        f"{accession_number}/{accession_number}.sra"
    ]
    
    for file_path in files_to_delete:
        if os.path.isfile(file_path):
            os.remove(file_path)
            log_message(f"Deleted {file_path}", level="info")
        else:
            log_message(f"File not found: {file_path}", level="warning")

def print_banner():
    # Generate ASCII banner using pyfiglet
    banner_text = pyfiglet.figlet_format("CANCER IMMUNOLOGY", font="slant")
    
    # Add a styled emoji banner
    emoji_banner = "🧬🔬🧪 CANCER RESEARCH - SNP PIPELINE 🧪🔬🧬"
    
    # Define the subtext with color highlights
    sub_text = "🔥 SNP Analysis Pipeline v3.0 🔥"

    # Print the enhanced banner
    print(colored(banner_text, "white"))
    print(colored(emoji_banner, "cyan", attrs=["bold"]))
    print(colored("=" * len(sub_text), "yellow"))
    print(colored(sub_text.center(len(sub_text)), "green", attrs=["bold"]))
    print(colored("=" * len(sub_text), "yellow"))

def run_command_with_progress(command, description="Processing", duration_estimate=100):
    with tqdm(total=duration_estimate, desc=description, unit="s", ncols=100) as pbar:
        process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        while process.poll() is None:
            time.sleep(1)  # Simulate work being done
            pbar.update(1)
        
        stdout, stderr = process.communicate()
        if process.returncode != 0:
            raise subprocess.CalledProcessError(process.returncode, command, output=stdout, stderr=stderr)
        
        pbar.update(duration_estimate - pbar.n)

def main():
    # Print banner
    print_banner()

    bwa_base_path = "/usr/local/bin/bwa/"
    bowtie_base_path = "/usr/local/bin/bowtie/"
    fastp_path = "/usr/bin/fastp"

    # Step 1: Detect and select the accession list file (.txt)
    accession_list_file = detect_accession_list_file()

    print()  # Line space after step 1

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

    incomplete_accessions = []
    for accession_number in accession_numbers:
        incomplete_chromosomes = []
        total_chromosomes = len(chromosomes_list)
        
        for chromosome in chromosomes_list:
            final_vcf_file = f"{accession_number}/{accession_number}_mapped_{chromosome}.var.-final.vcf"
            if not os.path.isfile(final_vcf_file) or os.path.getsize(final_vcf_file) == 0:
                incomplete_chromosomes.append(chromosome)
        
        if incomplete_chromosomes:
            incomplete_accessions.append({
                'accession': accession_number,
                'incomplete_count': len(incomplete_chromosomes),
                'total_count': total_chromosomes
            })

    total_incomplete = len(incomplete_accessions)
    log_message(f"Total incomplete accession numbers found: {total_incomplete}", level="info")

    for accession in incomplete_accessions:
        log_message(f"Accession {accession['accession']} has {accession['incomplete_count']} incomplete chromosomes out of {accession['total_count']}.", level="info")

    print()  # Line space after showing incomplete details

    if total_incomplete > 0:
        print(f"3. How many accession numbers would you like to analyze? (1-{total_incomplete}):", end=" ")
        num_to_analyze = int(input().strip())
        accession_numbers_to_analyze = incomplete_accessions[:num_to_analyze]
    else:
        log_message("All accession numbers are complete. Nothing to analyze.", level="success")
        return

    print()  # Line space after choosing how many accessions to analyze

    # Variable to track if trimmed files were already found for an accession number
    trimmed_files_found = {}

    for chromosome in chromosomes_list:
        for accession_info in accession_numbers_to_analyze:
            accession_number = accession_info['accession']
            ensure_directory(accession_number)

            # Detect if trimmed files already exist for this accession number
            if accession_number not in trimmed_files_found:
                trimmed_file_single = f"{accession_number}/{accession_number}_trimmed.fq"
                trimmed_file_paired_1 = f"{accession_number}/{accession_number}_1_trimmed.fq"
                trimmed_file_paired_2 = f"{accession_number}/{accession_number}_2_trimmed.fq"

                if os.path.isfile(trimmed_file_single) or (os.path.isfile(trimmed_file_paired_1) and os.path.isfile(trimmed_file_paired_2)):
                    log_message(f"Trimmed files already exist for {accession_number}. Skipping download and trimming.", level="warning")
                    read_type = detect_read_type(accession_number)
                    trimmed_files_found[accession_number] = True
                else:
                    # If no trimmed files are found, download and trim the FASTQ files
                    prefetch_and_convert(accession_number, threads)
                    read_type = detect_read_type(accession_number)
                    trim_reads(accession_number, read_type, fastp_path)
            else:
                # If trimmed files were found previously, skip download and trimming
                log_message(f"Trimmed files already detected earlier for {accession_number}. Skipping FASTQ check.", level="info")
                read_type = detect_read_type(accession_number)

            # Line space after prefetching or skipping trimming
            print()

            final_vcf_file = f"{accession_number}/{accession_number}_mapped_{chromosome}.var.-final.vcf"
            if os.path.isfile(final_vcf_file) and os.path.getsize(final_vcf_file) > 0:
                log_message(f"Final VCF file for chromosome {chromosome} already exists. Skipping analysis.", level="info")
                continue

            log_message(f"Processing chromosome {chromosome} for {accession_number}...", level="info")

            start_time = time.time()

            bwa_chrom_path = f"{bwa_base_path}{chromosome}_bwa_ind/Homo_sapiens.GRCh38.dna.chromosome.{chromosome}.fa" if chromosome != 'hg38' else f"{bwa_base_path}hg38/GRCh38_reference.fa"
            bowtie_index_path = f"{bowtie_base_path}{chromosome}_bowtie_ind/bowtie" if chromosome != 'hg38' else f"{bowtie_base_path}hg38/bowtie"

            input_file_for_mapping = trimmed_file_single if read_type == '1' else trimmed_file_paired_1

            log_message(f"\nMapping {accession_number} reads using Bowtie2 for chromosome {chromosome}...", level="info")
            
            bowtie_command = f"bowtie2 --very-fast-local -x {bowtie_index_path} -U {input_file_for_mapping} -S {accession_number}/{accession_number}_mapped_{chromosome}.sam"
            run_command_with_progress(bowtie_command, description=f"Mapping {accession_number} for chr {chromosome}", duration_estimate=120)

            log_message(f"Converting SAM to BAM for chromosome {chromosome}...", level="info")
            run_command(f"samtools view -S -b {accession_number}/{accession_number}_mapped_{chromosome}.sam > {accession_number}/{accession_number}_mapped_{chromosome}.bam")

            log_message(f"Sorting BAM for chromosome {chromosome}...", level="blue")
            run_command(f"samtools sort {accession_number}/{accession_number}_mapped_{chromosome}.bam > {accession_number}/{accession_number}_mapped_{chromosome}.sorted.bam")

            log_message(f"Generating VCF for chromosome {chromosome}...", level="blue")
            run_command(f"bcftools mpileup -f {bwa_chrom_path} {accession_number}/{accession_number}_mapped_{chromosome}.sorted.bam | bcftools call -mv -Ob -o {accession_number}/{accession_number}_mapped_{chromosome}.raw.bcf")

            log_message(f"Finalizing VCF for chromosome {chromosome}...", level="magenta")
            run_command(f"bcftools view {accession_number}/{accession_number}_mapped_{chromosome}.raw.bcf | vcfutils.pl varFilter - > {final_vcf_file}")

            delete_intermediate_files(accession_number, chromosome)

            end_time = time.time()
            elapsed_time = end_time - start_time
            log_message(f"Time elapsed for chromosome {chromosome}: {elapsed_time:.2f} seconds", level="info")

            print()

if __name__ == "__main__":
    main()