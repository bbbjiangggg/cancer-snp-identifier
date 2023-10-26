
import os
import shutil
import subprocess
from pathlib import Path
import signal
import importlib

def process_sra_sequence(sra, trim_path, truseq3_path, bowtie_index_path, bwa_chrom_path):
    trimmed_file = f"SRR_one/SRR_one_trimmed.fq.gz"
    if not os.path.isfile(trimmed_file):
        print("\nDownloading number sequence SRR_one from SRA... ")
        subprocess.run(["fastq-dump", sra])

        if os.path.isdir("SRR_one"):
            shutil.rmtree("SRR_one")
        os.makedirs("SRR_one")
        shutil.move(f"{sra}.fastq", "SRR_one")

        print("\nRunning fastqc on SRR_one... ")
        subprocess.run(["fastqc", f"SRR_one/{sra}.fastq"])

        print("\nTrimming SRR_one... ")
        subprocess.run(["java", "-jar", trim_path, "SE", f"SRR_one/{sra}.fastq", "SRR_one/SRR_one_trimmed.fq.gz", 
                        f"ILLUMINACLIP:{truseq3_path}:2:30:10", "SLIDINGWINDOW:4:20", "MINLEN:35"])

        print("\nRunning fastqc on trimmed SRR_one... ")
        subprocess.run(["fastqc", "SRR_one/SRR_one_trimmed.fq.gz"])
    else:
        print("\nTrimmed file already exists. Skipping download, trimming, and quality check...")

    print("\nMapping SRR_one reads using Bowtie2... ")
    subprocess.run(["bowtie2", "--very-fast-local", "-x", bowtie_index_path, "SRR_one/SRR_one_trimmed.fq.gz", 
                    "-S", "SRR_one/SRR_one_mapped.sam"])

    subprocess.run(["samtools", "view", "-S", "-b", "SRR_one/SRR_one_mapped.sam"], stdout=open("SRR_one/SRR_one_mapped.bam", "w"))

    print("\nSorting using Samtools... ")
    subprocess.run(["samtools", "sort", "SRR_one/SRR_one_mapped.bam"], stdout=open("SRR_one/SRR_one_mapped.sorted.bam", "w"))

    print("\nSummarizing the base calls (mpileup)... ")
    subprocess.run(["bcftools", "mpileup", "-f", bwa_chrom_path, "SRR_one/SRR_one_mapped.sorted.bam"], stdout=subprocess.PIPE)
    subprocess.run(["bcftools", "call", "-mv", "-Ob", "-o", "SRR_one/SRR_one_mapped.raw.bcf"])

    print("\nFinalizing VCF... ")
    subprocess.run(["bcftools", "view", "SRR_one/SRR_one_mapped.raw.bcf"], stdout=subprocess.PIPE)
    subprocess.run(["vcfutils.pl", "varFilter", "-"], stdout=open("SRR_one/SRR_one_mapped.var.-final.vcf", "w"))

    # Cleaning up
    for file_to_delete in ["SRR_one/SRR_one.fastq", "SRR_one/SRR_one_mapped.sam", "SRR_one/SRR_one_mapped.bam", 
                           "SRR_one/SRR_one_mapped.sorted.bam", "SRR_one/SRR_one_mapped.raw.bcf", "SRR_one/SRR_one_fastqc.zip", 
                           "SRR_one/SRR_one_trimmed_fastqc.zip"]:
        if os.path.exists(file_to_delete):
            os.remove(file_to_delete)


# Additional content from your script...

# ...

# Integration with existing script
for index, sra in enumerate(sra_list):
    print(f"\nProcessing {sra}...")
    process_sra_sequence(sra, trim_path, truseq3_path, bowtie_index_path, bwa_chrom_path)

# ...
