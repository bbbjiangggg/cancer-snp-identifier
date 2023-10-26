import os
import subprocess

def process_accession(accession, chromosomes_list):
    # Define the base paths for BWA, Bowtie, and Trimmomatic
    bwa_base_path = "/usr/local/bin/bwa/"
    bowtie_base_path = "/usr/local/bin/bowtie/"
    trimmomatic_path = "/usr/local/bin/Trimmomatic-0.39/trimmomatic-0.39.jar"
    truseq3_path = "/usr/local/bin/Trimmomatic-0.39/adapters/TruSeq3-SE.fa"

    # Define the paths of the potential trimmed file and the original fastq file
    trimmed_file = f"{accession}/{accession}_trimmed.fq.gz"
    fastq_file = f"{accession}/{accession}.fastq"

    # Check if the trimmed file already exists
    if not os.path.isfile(trimmed_file):
        print(f"\n\033[1;35mDownloading sequence {accession} from SRA...\033[0m")
        subprocess.run(["fastq-dump", accession])

        if os.path.isdir(accession):
            subprocess.run(["rm", "-r", accession])

        os.makedirs(accession, exist_ok=True)
        subprocess.run(["mv", f"{accession}.fastq", accession])

        print(f"\n\033[1;35mRunning fastqc on {accession}...\033[0m")
        subprocess.run(["fastqc", fastq_file])

        print(f"\n\033[1;35mTrimming {accession}...\033[0m")
        subprocess.run(["java", "-jar", trimmomatic_path, "SE", "-phred33", fastq_file, trimmed_file, "ILLUMINACLIP:" + truseq3_path + ":2:30:10", "SLIDINGWINDOW:4:20", "MINLEN:35"])

        print(f"\n\033[1;35mRunning fastqc on trimmed {accession}...\033[0m")
        subprocess.run(["fastqc", trimmed_file])
    else:
        print(f"\n\033[1;32mTrimmed file for {accession} already exists. Skipping download, trimming, and quality check...\033[0m")

    # Loop over each chromosome and perform the processing
    for chromosome in chromosomes_list:
        # Construct the paths
        bwa_chrom_path = f"{bwa_base_path}{chromosome}_bwa_ind/Homo_sapiens.GRCh38.dna.chromosome.{chromosome}.fa"
        bowtie_index_path = f"{bowtie_base_path}{chromosome}_bowtie_ind/bowtie"

        # Print the paths
        print(f"\nPaths for chromosome {chromosome}:")
        print("BWA Chromosome Path:", bwa_chrom_path)
        print("Bowtie Index Path:", bowtie_index_path)

        print(f"\n\033[1;35mMapping {accession} reads using Bowtie2 for chromosome {chromosome}...\033[0m")
        subprocess.run(["bowtie2", "--very-fast-local", "-x", bowtie_index_path, "-U", trimmed_file, "-S", f"{accession}/{accession}_{chromosome}_mapped.sam"])

        subprocess.run(["samtools", "view", "-S", "-b", f"{accession}/{accession}_{chromosome}_mapped.sam", "-o", f"{accession}/{accession}_{chromosome}_mapped.bam"])
        subprocess.run(["samtools", "sort", f"{accession}/{accession}_{chromosome}_mapped.bam", "-o", f"{accession}/{accession}_{chromosome}_mapped.sorted.bam"])

        print(f"\n\033[1;35mIndexing BAM file for chromosome {chromosome}...\033[0m")
        subprocess.run(["samtools", "index", f"{accession}/{accession}_{chromosome}_mapped.sorted.bam"])

        print(f"\n\033[1;35mSummarizing the base calls (mpileup) for chromosome {chromosome}...\033[0m")
        subprocess.run(["bcftools", "mpileup", "-f", bwa_chrom_path, f"{accession}/{accession}_{chromosome}_mapped.sorted.bam", "-o", f"{accession}/{accession}_{chromosome}_mapped.raw.bcf"])
        subprocess.run(["bcftools", "call", "-mv", "-Ob", "-o", f"{accession}/{accession}_{chromosome}_mapped.raw.bcf", f"{accession}/{accession}_{chromosome}_mapped.raw.bcf"])

        print(f"\n\033[1;35mFinalizing VCF for chromosome {chromosome}...\033[0m")
        subprocess.run(["bcftools", "view", f"{accession}/{accession}_{chromosome}_mapped.raw.bcf", "-o", f"{accession}/{accession}_{chromosome}_mapped.var.-final.vcf"])

        # Remove intermediate files for this chromosome
        '''subprocess.run(["rm", "-f", f"{accession}/{accession}_{chromosome}_mapped.sam", f"{accession}/{accession}_{chromosome}_mapped.bam", f"{accession}/{accession}_{chromosome}_mapped.sorted.bam", f"{accession}/{accession}_{chromosome}_mapped.raw.bcf"])'''

if __name__ == "__main__":
    # Prompt the user to enter the chromosomes
    chromosomes_input = input("Please enter the chromosomes to be analyzed, separated by a comma: ")
    chromosomes_list = [chromosome.strip() for chromosome in chromosomes_input.split(',')]
    print("List of chromosomes to be analyzed:", chromosomes_list)

    # Prompt the user for the name of the text file containing accession numbers
    accession_file = input("\nPlease enter the name of the .txt file that contains the list of accession numbers: ")

    try:
        # Read the accession numbers from the file
        with open(accession_file, 'r') as f:
            accessions_list = [line.strip() for line in f if line.strip()]
    except Exception as e:
        print(f"An error occurred while reading the file: {e}")
        exit(1)

    # Prompt the user for the number of accession numbers to analyze
    num_accessions = int(input(f"\nHow many accession numbers from the list do you want to analyze (max {len(accessions_list)}): "))
    if num_accessions > len(accessions_list):
        print("The specified number is greater than the number of accession numbers in the list. Exiting.")
        exit(1)

    # Process the specified number of accession numbers
    for accession in accessions_list[:num_accessions]:
        print(f"\nProcessing {accession}...")
        process_accession(accession, chromosomes_list)
