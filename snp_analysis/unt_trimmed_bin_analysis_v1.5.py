import os
import subprocess
import sys
import pyfiglet

def run_command(command):
    try:
        subprocess.run(command, check=True, shell=True)
    except subprocess.CalledProcessError as e:
        print(f"An error occurred: {e}", file=sys.stderr)
        sys.exit(1)
    except KeyboardInterrupt:
        print("\nAnalysis interrupted by user. Exiting.")
        sys.exit(1)

def print_chromosome_paths(chromosomes_list, bwa_base_path, bowtie_base_path):
    for chromosome in chromosomes_list:
        bwa_chrom_path = f"{bwa_base_path}{chromosome}_bwa_ind/Homo_sapiens.GRCh38.dna.chromosome.{chromosome}.fa"
        bowtie_index_path = f"{bowtie_base_path}{chromosome}_bowtie_ind/bowtie"
        print(f"\nPaths for chromosome {chromosome}:")
        print("BWA Chromosome Path:", bwa_chrom_path)
        print("Bowtie Index Path:", bowtie_index_path)

def read_accession_numbers(file_path):
    try:
        with open(file_path, 'r') as f:
            accession_numbers = [line.strip() for line in f if line.strip()]
        return accession_numbers
    except FileNotFoundError:
        print("The specified file was not found. Please check the file path and try again.")
        sys.exit(1)
    except Exception as e:
        print(f"An error occurred while reading the file: {e}")
        sys.exit(1)

def is_file_empty(file_path):
    return os.path.isfile(file_path) and os.path.getsize(file_path) == 0

def delete_intermediate_files(accession_number, chromosome):
    intermediate_files = [
        f"{accession_number}/{accession_number}.fastq",
        f"{accession_number}/{accession_number}_mapped_{chromosome}.sam",
        f"{accession_number}/{accession_number}_mapped_{chromosome}.bam",
        f"{accession_number}/{accession_number}_mapped_{chromosome}.raw.bcf",
        f"{accession_number}/{accession_number}_fastqc.zip",
        f"{accession_number}/{accession_number}_trimmed_fastqc.zip"
    ]
    for file_path in intermediate_files:
        if os.path.isfile(file_path):
            os.remove(file_path)
            print(f"Deleted {file_path}")

def main():
    text = "CANCER IMMUNOLOGY"
    font = "banner3-D"
    terminal_width = os.get_terminal_size().columns
    f = pyfiglet.Figlet(font=font, width=terminal_width)
    logo = f.renderText(text)
    print(logo.center(terminal_width))

    bwa_base_path = "/usr/local/bin/bwa/"
    bowtie_base_path = "/usr/local/bin/bowtie/"
    trimmomatic_path = "/usr/local/bin/Trimmomatic-0.39/trimmomatic-0.39.jar"
    truseq3_path = "/usr/local/bin/Trimmomatic-0.39/adapters/TruSeq3-SE.fa"

    user_email = input("Please enter your email address to receive a notification once the analysis is complete: ").strip()
    job_title = input("Please enter a job title for this analysis: ").strip()
    accession_list_file = input("Please enter the path to the accession list file: ").strip()
    
    accession_numbers = read_accession_numbers(accession_list_file)
    
    print(f"\nTotal accession numbers found: {len(accession_numbers)}")
    num_to_analyze = int(input("How many accession numbers do you want to analyze? "))
    accession_numbers_to_analyze = accession_numbers[:num_to_analyze]
    
    all_chromosomes = [str(i) for i in range(1, 23)] + ['X', 'Y', 'hg38']
    chromosomes_input = input("Please enter the chromosomes to be analyzed, separated by a comma, or type 'all' to analyze all chromosomes: ")
    if chromosomes_input.lower() == 'all':
        chromosomes_list = all_chromosomes
    else:
        chromosomes_list = [chromosome.strip() for chromosome in chromosomes_input.split(',')]
    
    print("List of chromosomes to be analyzed:", chromosomes_list)
    print_chromosome_paths(chromosomes_list, bwa_base_path, bowtie_base_path)

    for accession_number in accession_numbers_to_analyze:
        trimmed_file = f"{accession_number}/{accession_number}_trimmed.fq.gz"
        if not os.path.isfile(trimmed_file):
            print(f"\n\033[1;35mDownloading number sequence {accession_number} from SRA...\033[0m ")
            run_command(f"fastq-dump {accession_number}")

            if os.path.isdir(accession_number):
                os.rmdir(accession_number)
            os.makedirs(accession_number, exist_ok=True)
            os.rename(f"{accession_number}.fastq", f"{accession_number}/{accession_number}.fastq")

            print(f"\n\033[1;35mRunning fastqc on {accession_number}...\033[0m ")
            run_command(f"fastqc {accession_number}/{accession_number}.fastq")

            print(f"\n\033[1;35mTrimming {accession_number}...\033[0m ")
            trim_command = f"java -jar {trimmomatic_path} SE -phred33 {accession_number}/{accession_number}.fastq {trimmed_file} ILLUMINACLIP:{truseq3_path}:2:30:10 SLIDINGWINDOW:4:20 MINLEN:35"
            run_command(trim_command)

            print(f"\n\033[1;35mRunning fastqc on trimmed {accession_number}...\033[0m ")
            run_command(f"fastqc {trimmed_file}")
        else:
            print("\n\033[1;32mTrimmed file already exists. Skipping download, trimming, and quality check...\033[0m")

        for chromosome in chromosomes_list:
            final_vcf_file = f"{accession_number}/{accession_number}_mapped_{chromosome}.var.-final.vcf"
            if os.path.isfile(final_vcf_file) and not is_file_empty(final_vcf_file):
                print(f"\n\033[1;32mVCF file for {accession_number}, chromosome {chromosome} already exists. Skipping analysis...\033[0m")
                continue
            elif is_file_empty(final_vcf_file):
                print(f"\n\033[1;33mVCF file for {accession_number}, chromosome {chromosome} is empty. Deleting and adding to analysis...\033[0m")
                os.remove(final_vcf_file)

            bwa_chrom_path = f"{bwa_base_path}{chromosome}_bwa_ind/Homo_sapiens.GRCh38.dna.chromosome.{chromosome}.fa"
            bowtie_index_path = f"{bowtie_base_path}{chromosome}_bowtie_ind/bowtie"

            print(f"\n\033[1;35mMapping {accession_number} reads using Bowtie2 for chromosome {chromosome}...\033[0m ")
            run_command(f"bowtie2 --very-fast-local -x {bowtie_index_path} {trimmed_file} -S {accession_number}/{accession_number}_mapped_{chromosome}.sam")

            run_command(f"samtools view -S -b {accession_number}/{accession_number}_mapped_{chromosome}.sam > {accession_number}/{accession_number}_mapped_{chromosome}.bam")

            print("\n\033[1;35mSorting using Samtools...\033[0m ")
            run_command(f"samtools sort {accession_number}/{accession_number}_mapped_{chromosome}.bam > {accession_number}/{accession_number}_mapped_{chromosome}.sorted.bam")

            print("\n\033[1;35mSummarizing the base calls (mpileup)...\033[0m ")
            run_command(f"bcftools mpileup -f {bwa_chrom_path} {accession_number}/{accession_number}_mapped_{chromosome}.sorted.bam | bcftools call -mv -Ob -o {accession_number}/{accession_number}_mapped_{chromosome}.raw.bcf")

            print("\n\033[1;35mFinalizing VCF...\033[0m ")
            run_command(f"bcftools view {accession_number}/{accession_number}_mapped_{chromosome}.raw.bcf | vcfutils.pl varFilter - > {final_vcf_file}")

            # Delete intermediate files to save disk space
            delete_intermediate_files(accession_number, chromosome)

    send_email_command = f'sendemail -f sudoroot1775@outlook.com -t {user_email} -u "{job_title}_Analysis Done" -m "Ready to receive information for the next analysis." -s smtp-mail.outlook.com:587 -o tls=yes -xu sudoroot1775@outlook.com -xp ydAEwVVu2s7uENC'
    os.system(send_email_command)
    print("Analysis complete! A notification email has been sent.")

if __name__ == "__main__":
    main()
