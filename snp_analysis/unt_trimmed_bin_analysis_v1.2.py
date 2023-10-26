import os
import subprocess
import sys

def run_command(command):
    try:
        subprocess.run(command, check=True, shell=True)
    except subprocess.CalledProcessError as e:
        print(f"An error occurred: {e}", file=sys.stderr)
        sys.exit(1)

def print_chromosome_paths(chromosomes_list, bwa_base_path, bowtie_base_path):
    # Print the paths for each chromosome
    for chromosome in chromosomes_list:
        # Construct the paths
        bwa_chrom_path = f"{bwa_base_path}{chromosome}_bwa_ind/Homo_sapiens.GRCh38.dna.chromosome.{chromosome}.fa"
        bowtie_index_path = f"{bowtie_base_path}{chromosome}_bowtie_ind/bowtie"
    
        # Print the paths
        print(f"\nPaths for chromosome {chromosome}:")
        print("BWA Chromosome Path:", bwa_chrom_path)
        print("Bowtie Index Path:", bowtie_index_path)

def main():
    # Define the base paths
    bwa_base_path = "/usr/local/bin/bwa/"
    bowtie_base_path = "/usr/local/bin/bowtie/"
    trimmomatic_path = "/usr/local/bin/Trimmomatic-0.39/trimmomatic-0.39.jar"
    truseq3_path = "/usr/local/bin/Trimmomatic-0.39/adapters/TruSeq3-SE.fa"

    # Prompt the user to enter the accession number and chromosomes
    accession_number = input("Please enter the accession number: ").strip()
    chromosomes_input = input("Please enter the chromosomes to be analyzed, separated by a comma: ")
    
    # Split the input string by comma and remove leading/trailing whitespace
    chromosomes_list = [chromosome.strip() for chromosome in chromosomes_input.split(',')]
    
    # Print the list of chromosomes and their paths
    print("List of chromosomes to be analyzed:", chromosomes_list)
    print_chromosome_paths(chromosomes_list, bwa_base_path, bowtie_base_path)
    
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

    print(f"\n\033[1;35mMapping {accession_number} reads using Bowtie2...\033[0m ")
    run_command(f"bowtie2 --very-fast-local -x bowtie_index_path {trimmed_file} -S {accession_number}/{accession_number}_mapped.sam")

    run_command(f"samtools view -S -b {accession_number}/{accession_number}_mapped.sam > {accession_number}/{accession_number}_mapped.bam")

    print("\n\033[1;35mSorting using Samtools...\033[0m ")
    run_command(f"samtools sort {accession_number}/{accession_number}_mapped.bam > {accession_number}/{accession_number}_mapped.sorted.bam")

    print("\n\033[1;35mSummarizing the base calls (mpileup)...\033[0m ")
    run_command(f"bcftools mpileup -f bwa_chrom_path {accession_number}/{accession_number}_mapped.sorted.bam | bcftools call -mv -Ob -o {accession_number}/{accession_number}_mapped.raw.bcf")

    print("\n\033[1;35mFinalizing VCF...\033[0m ")
    run_command(f"bcftools view {accession_number}/{accession_number}_mapped.raw.bcf | vcfutils.pl varFilter - > {accession_number}/{accession_number}_mapped.var.-final.vcf")

if __name__ == "__main__":
    main()
