import os 
import subprocess
import sys
import shutil
import pyfiglet
from termcolor import colored
import time  # Importing time module

# Function to check and install missing packages
def check_and_install(package):
    try:
        __import__(package)
    except ImportError:
        print(f"Package {package} not found. Installing...")
        subprocess.check_call([sys.executable, "-m", "pip", "install", package])

# Check and install required packages
required_packages = ['os', 'subprocess', 'sys', 'shutil', 'pyfiglet', 'termcolor', 'time']
for package in required_packages:
    check_and_install(package)

# Declare these variables as global so they can be accessed in functions
global user_email, job_title

def run_command(command):
    try:
        subprocess.run(command, check=True, shell=True)
    except subprocess.CalledProcessError as e:
        log_message(f"An error occurred: {e}", level="error")
        sys.exit(1)
    except KeyboardInterrupt:
        log_message("Analysis interrupted by user. Exiting.", level="error")
        sys.exit(1)

def ensure_directory(path):
    # Ensure the directory exists or create it
    if not os.path.exists(path):
        os.makedirs(path, exist_ok=True)
        log_message(f"Directory {path} created.", level="success")
    else:
        log_message(f"Directory {path} already exists.", level="warning")

def print_chromosome_paths(chromosomes_list, bwa_base_path, bowtie_base_path):
    for chromosome in chromosomes_list:
        bwa_chrom_path = f"{bwa_base_path}{chromosome}_bwa_ind/Homo_sapiens.GRCh38.dna.chromosome.{chromosome}.fa" if chromosome != 'hg38' else f"{bwa_base_path}hg38/GRCh38_reference.fa"
        bowtie_index_path = f"{bowtie_base_path}{chromosome}_bowtie_ind/bowtie" if chromosome != 'hg38' else f"{bowtie_base_path}hg38/bowtie"
        log_message(f"Paths for chromosome {chromosome}:", level="info")
        log_message("BWA Chromosome Path: ", normal_part=bwa_chrom_path, level="info")
        log_message("Bowtie Index Path: ", normal_part=bowtie_index_path, level="info")

def read_accession_numbers(file_path):
    try:
        with open(file_path, 'r') as f:
            accession_numbers = [line.strip() for line in f if line.strip()]
        return accession_numbers
    except FileNotFoundError:
        log_message("The specified file was not found. Please check the file path and try again.", level="error")
        sys.exit(1)
    except Exception as e:
        log_message(f"An error occurred while reading the file: {e}", level="error")
        sys.exit(1)

def is_file_empty(file_path):
    return os.path.isfile(file_path) and os.path.getsize(file_path) == 0

def delete_intermediate_files(accession_number, chromosome):
    intermediate_files = [
        f"{accession_number}/{accession_number}.fastq",
        f"{accession_number}/{accession_number}_mapped_{chromosome}.sam",
        f"{accession_number}/{accession_number}_mapped_{chromosome}.bam",
        f"{accession_number}/{accession_number}_mapped_{chromosome}.raw.bcf",
        f"{accession_number}/{accession_number}_fastqc.zip"
    ]
    for file_path in intermediate_files:
        if os.path.isfile(file_path):
            os.remove(file_path)
            log_message(f"Deleted {file_path}", level="success")

def clear_directory(directory):
    if os.path.isdir(directory):
        shutil.rmtree(directory)
        os.makedirs(directory, exist_ok=True)

def get_verified_path(prompt_message):
    while True:
        path = input(colored(prompt_message, "cyan")).strip()
        if os.path.exists(path):
            return path
        else:
            log_message("The provided path does not exist. Please try again.", level="warning")

def delete_existing_fastq_files(accession_number, read_type):
    if read_type == '1':
        fastq_file = f"{accession_number}/{accession_number}.fastq"
        if os.path.isfile(fastq_file):
            os.remove(fastq_file)
            log_message(f"Deleted existing file {fastq_file}", level="warning")
    elif read_type == '2':
        fastq_file_1 = f"{accession_number}/{accession_number}_1.fastq"
        fastq_file_2 = f"{accession_number}/{accession_number}_2.fastq"
        if os.path.isfile(fastq_file_1):
            os.remove(fastq_file_1)
            log_message(f"Deleted existing file {fastq_file_1}", level="warning")
        if os.path.isfile(fastq_file_2):
            os.remove(fastq_file_2)
            log_message(f"Deleted existing file {fastq_file_2}", level="warning")

def log_message(message, level="normal", normal_part=""):
    color_map = {
        "info": "\033[95m",       # Magenta
        "warning": "\033[93m",    # Yellow
        "error": "\033[91m",      # Red
        "success": "\033[92m",    # Green
        "normal": "\033[0m",      # Reset to default
        "blue": "\033[94m"        # Blue
    }
    color = color_map.get(level, "\033[0m")
    if normal_part:
        print(f"{color}{message}\033[0m{normal_part}")
    else:
        print(f"{color}{message}\033[0m")

# Function to automatically detect the .txt file and ask the user for confirmation
def detect_accession_list_file():
    txt_files = [f for f in os.listdir() if f.endswith('.txt')]
    if txt_files:
        accession_list_file = txt_files[0]
        log_message(f"Detected accession list file: {accession_list_file}. Is this correct? (yes/y or no/n)", level="info")
        user_input = input().strip().lower()
        if user_input in ['yes', 'y']:
            return accession_list_file
        else:
            return get_verified_path("Please enter the path to the accession list file: ")
    else:
        return get_verified_path("No .txt file detected. Please enter the path to the accession list file: ")

# Function to check if a file is compressed
def is_compressed(file_path):
    return file_path.endswith('.gz')

def main():
    global user_email, job_title  # Use global declaration
    text = "CANCER IMMUNOLOGY v2.0 - SNP Analysis Pipeline"
    terminal_width = os.get_terminal_size().columns
    print("=" * terminal_width)
    print(text.center(terminal_width))
    print("=" * terminal_width)

    # Define base paths for BWA and Bowtie2
    bwa_base_path = "/usr/local/bin/bwa/"
    bowtie_base_path = "/usr/local/bin/bowtie/"
    fastp_path = "/usr/bin/fastp"

    # Verify paths for BWA and Bowtie2 indexes
    if not os.path.exists(bwa_base_path):
        bwa_base_path = get_verified_path("1. BWA base path not found. Please enter the correct BWA base path: ")
    if not os.path.exists(bowtie_base_path):
        bowtie_base_path = get_verified_path("2. Bowtie base path not found. Please enter the correct Bowtie base path: ")

    # Step 1: Detect and verify the accession list file
    accession_list_file = detect_accession_list_file()

    # Step 2: Are the reads single-end or paired-end?
    print("2. Are the reads single-end (1) or paired-end (2)? Enter 1 or 2:", end=" ")
    read_type = input().strip()
    print()  # Line space after input

    # Step 3: Enter the number of threads to use
    available_threads = os.cpu_count()
    recommended_threads = max(1, available_threads // 2)
    print(f"3. Enter the number of threads to use (Available: {available_threads}, Recommended: {recommended_threads}):", end=" ")
    threads = input().strip()
    if not threads.isdigit() or int(threads) < 1 or int(threads) > available_threads:
        log_message(f"Invalid input. Using recommended number of threads: {recommended_threads}.", level="warning")
        threads = str(recommended_threads)
    print()  # Line space after input

    # Step 4: Enter chromosomes to analyze
    print("4. Enter chromosomes to analyze (comma-separated) or 'all' for all chromosomes:", end=" ")
    chromosomes_input = input().strip()
    print()  # Line space after input

    accession_numbers = read_accession_numbers(accession_list_file)
    total_accession_numbers = len(accession_numbers)
    all_chromosomes = [str(i) for i in range(1, 23)] + ['X', 'Y']

    # Display total accession numbers found
    log_message(f"Total accession numbers found: {total_accession_numbers}", level="info")
    print()  # Line space after input

    chromosomes_list = all_chromosomes if chromosomes_input.lower() == 'all' else [chromosome.strip() for chromosome in chromosomes_input.split(',')]

    accession_numbers_to_analyze = []

    for accession_number in accession_numbers:
        analysis_needed = False
        for chromosome in chromosomes_list:
            final_vcf_file = f"{accession_number}/{accession_number}_mapped_{chromosome}.var.-final.vcf"
            if not os.path.isfile(final_vcf_file) or is_file_empty(final_vcf_file):
                analysis_needed = True
                log_message(f"Final VCF file for {accession_number}, chromosome {chromosome} is missing or empty.", level="info")
                break
            else:
                log_message(f"VCF file for {accession_number}, chromosome {chromosome} already exists. Skipping analysis for this chromosome...", level="warning")

        if analysis_needed:
            accession_numbers_to_analyze.append(accession_number)

    # Step 5: How many accession numbers to analyze?
    print(f"5. How many accession numbers would you like to analyze? (1-{total_accession_numbers}):", end=" ")
    num_to_analyze = int(input().strip())
    print()  # Line space after input

    accession_numbers_to_analyze = accession_numbers_to_analyze[:num_to_analyze]

    log_message("List of chromosomes to be analyzed:", level="info", normal_part=str(chromosomes_list))
    print_chromosome_paths(chromosomes_list, bwa_base_path, bowtie_base_path)

    # Prepare data for all accession_numbers
    for accession_number in accession_numbers_to_analyze:
        # Ensure the output directory exists, or create it
        ensure_directory(accession_number)

        # Check if trimmed files exist before attempting any downloads or trimming
        trimmed_file_single = f"{accession_number}/{accession_number}_trimmed.fq"
        trimmed_file_single_gz = trimmed_file_single + ".gz"
        trimmed_file_paired_1 = f"{accession_number}/{accession_number}_1_trimmed.fq"
        trimmed_file_paired_1_gz = trimmed_file_paired_1 + ".gz"
        trimmed_file_paired_2 = f"{accession_number}/{accession_number}_2_trimmed.fq"
        trimmed_file_paired_2_gz = trimmed_file_paired_2 + ".gz"

        if (read_type == '1' and (os.path.isfile(trimmed_file_single) or os.path.isfile(trimmed_file_single_gz))) or \
           (read_type == '2' and ((os.path.isfile(trimmed_file_paired_1) and os.path.isfile(trimmed_file_paired_2)) or \
                                  (os.path.isfile(trimmed_file_paired_1_gz) and os.path.isfile(trimmed_file_paired_2_gz)))):
            log_message(f"Trimmed files for {accession_number} found. Skipping download and trimming.", level="warning")
        else:
            # Attempt to download the data using prefetch and fasterq-dump
            try:
                log_message(f"\nDownloading sequence number {accession_number} with prefetch...", level="info")
                run_command(f"prefetch {accession_number}")

                if read_type == '1':
                    log_message(f"\nConverting single-end sequence number {accession_number} with fasterq-dump...", level="info")
                    run_command(f"fasterq-dump {accession_number} --threads {threads} --progress")
                    # Move the resulting file to the corresponding directory
                    shutil.move(f"{accession_number}.fastq", f"{accession_number}/{accession_number}.fastq")
                elif read_type == '2':
                    log_message(f"\nConverting paired-end sequence number {accession_number} with fasterq-dump...", level="info")
                    run_command(f"fasterq-dump {accession_number} --split-files --threads {threads} --progress")
                    # Move the resulting files to the corresponding directory
                    shutil.move(f"{accession_number}_1.fastq", f"{accession_number}/{accession_number}_1.fastq")
                    shutil.move(f"{accession_number}_2.fastq", f"{accession_number}/{accession_number}_2.fastq")
                else:
                    log_message(f"Invalid read type. Skipping {accession_number}.", level="error")
                    continue

                log_message(f"FASTQ files for {accession_number} have been moved to the directory {accession_number}.", level="success")

                # Now, continue with trimming the files
                if read_type == '1':
                    log_message(f"\nTrimming {accession_number} with fastp...", level="info")
                    trim_command = f"{fastp_path} -i {accession_number}/{accession_number}.fastq -o {trimmed_file_single} --thread=4"
                    run_command(trim_command)
                elif read_type == '2':
                    log_message(f"Trimming {accession_number} with fastp...", level="info")
                    trim_command = f"{fastp_path} -i {accession_number}/{accession_number}_1.fastq -I {accession_number}/{accession_number}_2.fastq -o {trimmed_file_paired_1} -O {trimmed_file_paired_2} --thread=4"
                    run_command(trim_command)

                shutil.move("fastp.html", f"{accession_number}/fastp.html")
                shutil.move("fastp.json", f"{accession_number}/fastp.json")

                # Optionally compress the trimmed files to save space
                compress = input(colored("Would you like to compress the trimmed files to save space? (yes/y or no/n): ", "cyan")).strip().lower()
                if compress in ['yes', 'y']:
                    if read_type == '1':
                        run_command(f"gzip {trimmed_file_single}")
                    elif read_type == '2':
                        run_command(f"gzip {trimmed_file_paired_1}")
                        run_command(f"gzip {trimmed_file_paired_2}")

            except subprocess.CalledProcessError as e:
                log_message(f"An error occurred while processing {accession_number} with prefetch or fasterq-dump. Please check the log for details.", level="error")
                sys.exit(1)

    # Now process per chromosome
    for chromosome in chromosomes_list:
        for accession_number in accession_numbers_to_analyze:
            # Start a timer for each chromosome analysis
            start_time = time.time()

            final_vcf_file = f"{accession_number}/{accession_number}_mapped_{chromosome}.var.-final.vcf"

            if os.path.isfile(final_vcf_file) and not is_file_empty(final_vcf_file):
                log_message(f"VCF file for {accession_number}, chromosome {chromosome} already exists. Skipping analysis...", level="warning")
                continue
            elif is_file_empty(final_vcf_file):
                log_message(f"VCF file for {accession_number}, chromosome {chromosome} is empty. Deleting and adding to analysis...", level="warning")
                os.remove(final_vcf_file)

            # Set paths for BWA and Bowtie2 indices dynamically based on the chromosome
            if chromosome != 'hg38':
                bowtie_index_path = f"{bowtie_base_path}{chromosome}_bowtie_ind/bowtie"
                bwa_chrom_path = f"{bwa_base_path}{chromosome}_bwa_ind/Homo_sapiens.GRCh38.dna.chromosome.{chromosome}.fa"
            else:
                bowtie_index_path = f"{bowtie_base_path}hg38/bowtie"
                bwa_chrom_path = f"{bwa_base_path}hg38/GRCh38_reference.fa"

            # Determine if the trimmed files are compressed
            if read_type == '1':
                trimmed_file = f"{accession_number}/{accession_number}_trimmed.fq"
                if os.path.isfile(trimmed_file):
                    pass  # Uncompressed file exists
                elif os.path.isfile(trimmed_file + ".gz"):
                    trimmed_file += ".gz"  # Use compressed file
                else:
                    log_message(f"Trimmed file for {accession_number} not found.", level="error")
                    continue
                bowtie_input_option = f"-U {trimmed_file}"
            elif read_type == '2':
                trimmed_file_1 = f"{accession_number}/{accession_number}_1_trimmed.fq"
                trimmed_file_2 = f"{accession_number}/{accession_number}_2_trimmed.fq"
                if os.path.isfile(trimmed_file_1) and os.path.isfile(trimmed_file_2):
                    pass  # Uncompressed files exist
                elif os.path.isfile(trimmed_file_1 + ".gz") and os.path.isfile(trimmed_file_2 + ".gz"):
                    trimmed_file_1 += ".gz"  # Use compressed files
                    trimmed_file_2 += ".gz"
                else:
                    log_message(f"Trimmed files for {accession_number} not found.", level="error")
                    continue
                bowtie_input_option = f"-1 {trimmed_file_1} -2 {trimmed_file_2}"
            else:
                log_message("Invalid read type specified. Exiting.", level="error")
                sys.exit(1)

            # Remove '--gzip' from the bowtie_command
            bowtie_command = f"bowtie2 --very-fast-local -x {bowtie_index_path} {bowtie_input_option} --threads {threads} -S {accession_number}/{accession_number}_mapped_{chromosome}.sam"

            log_message(f"\nMapping {accession_number} reads using Bowtie2 for chromosome {chromosome}...", level="info")
            run_command(bowtie_command)

            run_command(f"samtools view -S -b {accession_number}/{accession_number}_mapped_{chromosome}.sam > {accession_number}/{accession_number}_mapped_{chromosome}.bam")

            log_message("\nSorting using Samtools...", level="blue")
            run_command(f"samtools sort {accession_number}/{accession_number}_mapped_{chromosome}.bam -o {accession_number}/{accession_number}_mapped_{chromosome}.sorted.bam")

            # Continue with the rest of your analysis pipeline...

            # Delete intermediate files to save disk space
            delete_intermediate_files(accession_number, chromosome)

            # End the timer for this chromosome and display elapsed time
            end_time = time.time()
            elapsed_time = end_time - start_time
            log_message(f"Time elapsed for accession number {accession_number}, chromosome {chromosome}: {elapsed_time:.2f} seconds", level="info")

if __name__ == "__main__":
    main()
