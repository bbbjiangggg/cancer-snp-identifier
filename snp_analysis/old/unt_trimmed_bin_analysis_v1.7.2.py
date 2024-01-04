import os
import subprocess
import sys
import shutil
import pyfiglet
import socket
import sendgrid
from sendgrid.helpers.mail import Mail, Email, To, Content
from termcolor import colored

def run_command(command):
    try:
        subprocess.run(command, check=True, shell=True)
    except subprocess.CalledProcessError as e:
        print(colored(f"An error occurred: {e}", "red"), file=sys.stderr)
        sys.exit(1)
    except KeyboardInterrupt:
        print(colored("\nAnalysis interrupted by user. Exiting.", "red"))
        sys.exit(1)

def print_chromosome_paths(chromosomes_list, bwa_base_path, bowtie_base_path):
    for chromosome in chromosomes_list:
        bwa_chrom_path = f"{bwa_base_path}{chromosome}_bwa_ind/Homo_sapiens.GRCh38.dna.chromosome.{chromosome}.fa" if chromosome != 'hg38' else f"{bwa_base_path}hg38/GRCh38_reference.fa"
        bowtie_index_path = f"{bowtie_base_path}{chromosome}_bowtie_ind/bowtie" if chromosome != 'hg38' else f"{bowtie_base_path}hg38/bowtie"
        print(colored(f"\nPaths for chromosome {chromosome}:", "cyan"))
        print(colored("BWA Chromosome Path:", "cyan"), bwa_chrom_path)
        print(colored("Bowtie Index Path:", "cyan"), bowtie_index_path)

def read_accession_numbers(file_path):
    try:
        with open(file_path, 'r') as f:
            accession_numbers = [line.strip() for line in f if line.strip()]
        return accession_numbers
    except FileNotFoundError:
        print(colored("The specified file was not found. Please check the file path and try again.", "red"))
        sys.exit(1)
    except Exception as e:
        print(colored(f"An error occurred while reading the file: {e}", "red"))
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
            print(colored(f"Deleted {file_path}", "green"))

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
            print(colored("The provided path does not exist. Please try again.", "yellow"))

def send_email_via_sendgrid(from_email, to_email, job_title, hostname, content):
    api_key = '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'  # Your SendGrid API Key
    sg = sendgrid.SendGridAPIClient(api_key=api_key)
    
    subject = f"{job_title} on {hostname} Analysis Completed"
    from_email = Email(from_email)
    to_email = To(to_email)
    content = Content("text/plain", content)

    mail = Mail(from_email, to_email, subject, content)

    response = sg.client.mail.send.post(request_body=mail.get())
    
    # Simplified user-friendly output based on the response
    if response.status_code == 202:
        print("Email sent successfully! Please check your inbox.")
    elif 400 <= response.status_code < 500:
        print("Failed to send email: There was a problem with the request. Please check the details and try again.")
    elif 500 <= response.status_code < 600:
        print("Failed to send email: Server error occurred. Please try again later.")
    else:
        print(f"Failed to send email: Unexpected status code received: {response.status_code}")

def main():
    text = "CANCER IMMUNOLOGY v1.7"
    font = "banner3-D"
    terminal_width = os.get_terminal_size().columns
    f = pyfiglet.Figlet(font=font, width=terminal_width)
    logo = f.renderText(text)
    print(logo.center(terminal_width))

    bwa_base_path = "/usr/local/bin/bwa/"
    bowtie_base_path = "/usr/local/bin/bowtie/"
    fastp_path = "/usr/bin/fastp"

    if not os.path.exists(bwa_base_path):
        bwa_base_path = get_verified_path("1. BWA base path not found. Please enter the correct BWA base path: ")
    if not os.path.exists(bowtie_base_path):
        bowtie_base_path = get_verified_path("2. Bowtie base path not found. Please enter the correct Bowtie base path: ")

    print(colored("1. Please enter your email address to receive a notification once the analysis is complete: ", "white"))
    user_email = input().strip()
    print(colored("2. Please enter a job title for this analysis: ", "white"))
    job_title = input().strip()
    print(colored("3. Please enter the path to the accession list file: ", "white"))
    accession_list_file = input().strip()

    print(colored("4. Are the reads single-end (1) or paired-end (2)? Enter 1 or 2: ", "white"))
    read_type = input().strip()

    accession_numbers = read_accession_numbers(accession_list_file)
    total_accession_numbers = len(accession_numbers)
    all_chromosomes = [str(i) for i in range(1, 23)] + ['X', 'Y']

    print(colored(f"\nTotal accession numbers in the file: {total_accession_numbers}", "cyan"))

    chromosomes_input = input(colored("5. Please enter the chromosomes to be analyzed, separated by a comma, or type 'all' to analyze all chromosomes: ", "white"))
    chromosomes_list = all_chromosomes if chromosomes_input.lower() == 'all' else [chromosome.strip() for chromosome in chromosomes_input.split(',')]

    accession_numbers_to_analyze = []

    for accession_number in accession_numbers:
        analysis_needed = False
        for chromosome in chromosomes_list:
            final_vcf_file = f"{accession_number}/{accession_number}_mapped_{chromosome}.var.-final.vcf"
            if not os.path.isfile(final_vcf_file) or is_file_empty(final_vcf_file):
                analysis_needed = True
                print(colored(f"Final VCF file for chromosome {chromosome} is missing.", "green"))
                break
            else:
                print(colored(f"VCF file for {accession_number}, chromosome {chromosome} already exists. Skipping analysis for this chromosome...", "yellow"))

        if analysis_needed:
            accession_numbers_to_analyze.append(accession_number)

    print(colored(f"Accession numbers pending analysis: {len(accession_numbers_to_analyze)}", "green"))

    if not accession_numbers_to_analyze:
        print(colored("No accession numbers left to analyze. Exiting.", "red"))
        sys.exit(0)

    num_to_analyze = int(input(colored("6. How many accession numbers do you want to analyze? ", "white")))
    accession_numbers_to_analyze = accession_numbers_to_analyze[:num_to_analyze]

    print(colored("List of chromosomes to be analyzed:", "cyan"), chromosomes_list)
    print_chromosome_paths(chromosomes_list, bwa_base_path, bowtie_base_path)

    for accession_number in accession_numbers_to_analyze:
        # Define file paths for single-end and paired-end reads
        trimmed_file_single = f"{accession_number}/{accession_number}_trimmed.fq.gz"
        trimmed_file_paired_1 = f"{accession_number}/{accession_number}_1_trimmed.fq.gz"
        trimmed_file_paired_2 = f"{accession_number}/{accession_number}_2_trimmed.fq.gz"

        # Download and trimming logic based on read type
        if read_type == '1' and not os.path.isfile(trimmed_file_single):
            print(colored(f"\nDownloading single-end sequence number {accession_number} from SRA...", "cyan"))
            run_command(f"fastq-dump {accession_number} -O {accession_number}/")
            print(colored(f"\nTrimming {accession_number} with fastp...", "cyan"))
            trim_command = f"{fastp_path} -i {accession_number}/{accession_number}.fastq -o {trimmed_file_single} --thread=4"
            run_command(trim_command)
            shutil.move("fastp.html", f"{accession_number}/fastp.html")
            shutil.move("fastp.json", f"{accession_number}/fastp.json")
        elif read_type == '2' and not (os.path.isfile(trimmed_file_paired_1) and os.path.isfile(trimmed_file_paired_2)):
            print(colored(f"\nDownloading paired-end sequence number {accession_number} from SRA...", "cyan"))
            run_command(f"fastq-dump --split-files {accession_number} -O {accession_number}/")
            print(colored(f"Trimming {accession_number} with fastp...", "cyan"))
            trim_command = f"{fastp_path} -i {accession_number}/{accession_number}_1.fastq -I {accession_number}/{accession_number}_2.fastq -o {trimmed_file_paired_1} -O {trimmed_file_paired_2} --thread=4"
            run_command(trim_command)
            shutil.move("fastp.html", f"{accession_number}/fastp.html")
            shutil.move("fastp.json", f"{accession_number}/fastp.json")
        else:
            print(colored(f"Trimmed files for {accession_number} found or invalid read type. Skipping download and trimming.", "yellow"))

        for chromosome in chromosomes_list:
            final_vcf_file = f"{accession_number}/{accession_number}_mapped_{chromosome}.var.-final.vcf"
            
            if os.path.isfile(final_vcf_file) and not is_file_empty(final_vcf_file):
                print(colored(f"VCF file for {accession_number}, chromosome {chromosome} already exists. Skipping analysis...", "yellow"))
                continue
            elif is_file_empty(final_vcf_file):
                print(colored(f"VCF file for {accession_number}, chromosome {chromosome} is empty. Deleting and adding to analysis...", "yellow"))
                os.remove(final_vcf_file)

            # Set paths for BWA and Bowtie2 indices
            bwa_chrom_path = f"{bwa_base_path}{chromosome}_bwa_ind/Homo_sapiens.GRCh38.dna.chromosome.{chromosome}.fa" if chromosome != 'hg38' else f"{bwa_base_path}hg38/GRCh38_reference.fa"
            bowtie_index_path = f"{bowtie_base_path}{chromosome}_bowtie_ind/bowtie" if chromosome != 'hg38' else f"{bowtie_base_path}hg38/bowtie"

            # Set the correct input file for mapping based on read type
            input_file_for_mapping = trimmed_file_single if read_type == '1' else trimmed_file_paired_1

            print(colored(f"\nMapping {accession_number} reads using Bowtie2 for chromosome {chromosome}...", "cyan"))
            bowtie_command = f"bowtie2 --very-fast-local -x {bowtie_index_path} -U {input_file_for_mapping} -S {accession_number}/{accession_number}_mapped_{chromosome}.sam"
            run_command(bowtie_command)

            run_command(f"samtools view -S -b {accession_number}/{accession_number}_mapped_{chromosome}.sam > {accession_number}/{accession_number}_mapped_{chromosome}.bam")

            print(colored("\nSorting using Samtools...", "magenta"))
            run_command(f"samtools sort {accession_number}/{accession_number}_mapped_{chromosome}.bam > {accession_number}/{accession_number}_mapped_{chromosome}.sorted.bam")

            print(colored("\nSummarizing the base calls (mpileup)...", "magenta"))
            run_command(f"bcftools mpileup -f {bwa_chrom_path} {accession_number}/{accession_number}_mapped_{chromosome}.sorted.bam | bcftools call -mv -Ob -o {accession_number}/{accession_number}_mapped_{chromosome}.raw.bcf")

            print(colored("\n\033[1;35mFinalizing VCF...\033[0m ", "magenta"))
            run_command(f"bcftools view {accession_number}/{accession_number}_mapped_{chromosome}.raw.bcf | vcfutils.pl varFilter - > {final_vcf_file}")

            # Delete intermediate files to save disk space
            delete_intermediate_files(accession_number, chromosome)

    # Email sending part replaced with SendGrid
    from_email = 'sudoroot1775@outlook.com'  # Your "from" email address
    to_email = user_email  # Use the previously input user_email
    job_title = job_title  # Use the previously input job_title
    hostname = socket.gethostname()  # Automatically detect the host computer name
    
    content = "The analysis has completed. Ready to receive information for the next analysis."

    send_email_via_sendgrid(from_email, to_email, job_title, hostname, content)

    

if __name__ == "__main__":
    main()
