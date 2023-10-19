#!/usr/bin/env python3

import os
import subprocess
from pathlib import Path

# Define color codes for consistent terminal output
RED = '\033[1;31m'
GREEN = '\033[1;32m'
MAGENTA = '\033[1;35m'
BLUE = '\033[1;34m'
RESET = '\033[0m'

# Function to execute shell commands and display colored output
def execute_command(command, message):
    print(message)
    process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()

    if process.returncode != 0:
        print(f"{RED}Error executing command: {command}{RESET}")
        print(stderr.decode())
    else:
        print(stdout.decode())

# Main script execution
if __name__ == "__main__":
    # Ask the user to input the chromosome number
    chrom_number = input("Please enter the chromosome number: ")

    # Construct the ref_chrom_path
    ref_chrom_path = f"/usr/local/bin/bwa/9_bwa_ind/Homo_sapiens.GRCh38.dna.chromosome.{chrom_number}.fa"

    # Check if the reference chromosome file exists
    if not Path(ref_chrom_path).is_file():
        print(f"{RED}The reference chromosome file does not exist at {ref_chrom_path}{RESET}")
        exit(1)

    # Scan for directories starting with ERR or SRR
    dirs = [d for d in os.listdir() if os.path.isdir(d) and (d.startswith("ERR") or d.startswith("SRR"))]

    for dir in dirs:
        print(f"{GREEN}Processing directory: {dir}{RESET}")

        # Command 0: Sort bam file
        cmd0 = f"samtools sort {dir}/{dir}_mapped.bam > {dir}/{dir}_mapped.sorted.bam"
        execute_command(cmd0, f"\n{MAGENTA}Sorting BAM file...{RESET}")

        # Command 1: Summarize base calls
        cmd1 = f"bcftools mpileup -f {ref_chrom_path} {dir}/{dir}_mapped.sorted.bam | bcftools call -mv -Ob -o {dir}/{dir}_mapped.raw.bcf"
        execute_command(cmd1, f"\n{MAGENTA}Summarizing the base calls (mpileup)...{RESET}")

        # Command 2: Finalize VCF
        cmd2 = f"bcftools view {dir}/{dir}_mapped.raw.bcf | vcfutils.pl varFilter - > {dir}/{dir}_mapped.var.-final.vcf"
        execute_command(cmd2, f"\n{MAGENTA}Finalizing VCF...{RESET}")

    print(f"{GREEN}Script execution completed.{RESET}")

