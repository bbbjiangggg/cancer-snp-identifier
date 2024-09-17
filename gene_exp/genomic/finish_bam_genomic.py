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
    # Scan for directories starting with ERR or SRR
    dirs = [d for d in os.listdir() if os.path.isdir(d) and (d.startswith("ERR") or d.startswith("SRR"))]

    for dir in dirs:
        print(f"{GREEN}Processing directory: {dir}{RESET}")

        # Replace 'ref_chrom_path' with the actual path to your reference chromosome file
        ref_chrom_path = "/usr/local/bin/bwa/hg38/GRCh38_reference.fa"

        # Command 1: Summarize base calls
        cmd1 = f"bcftools mpileup -f {ref_chrom_path} {dir}/{dir}_mapped_hg38.sorted.bam | bcftools call -mv -Ob -o {dir}/{dir}_mapped.raw.bcf"
        execute_command(cmd1, f"\n{MAGENTA}Summarizing the base calls (mpileup)...{RESET}")

        # Command 2: Finalize VCF
        cmd2 = f"bcftools view {dir}/{dir}_mapped.raw.bcf | vcfutils.pl varFilter - > {dir}/{dir}_mapped.var.-final.vcf"
        execute_command(cmd2, f"\n{MAGENTA}Finalizing VCF...{RESET}")

    print(f"{GREEN}Script execution completed.{RESET}")

