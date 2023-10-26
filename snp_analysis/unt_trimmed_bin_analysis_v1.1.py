#!/usr/bin/env python3

import os
import shutil
import subprocess
from pathlib import Path
import signal
import importlib

packages = ['os', 'shutil', 'subprocess', 'pathlib', 'signal']

# Check if packages are installed, install them if necessary
for package in packages:
    try:
        importlib.import_module(package)
    except ImportError:
        print(f'{package} is not installed. Installing...')
        subprocess.run(['pip', 'install', package])

# Define color codes
RED = '\033[1;31m'
GREEN = '\033[1;32m'
MAGENTA = '\033[1;35m'
BLUE = '\033[1;34m'
RESET = '\033[0m'

def replace_file_on_interrupt(sig, frame):
    # Remove the untrimmed_bash_sra_v1.2.txt file
    os.remove('untrimmed_bash_sra_v1.2.txt')
    print('\nRemoving untrimmed_bash_sra_v1.2.txt...')
    exit(1)

# Register the signal handler
signal.signal(signal.SIGINT, replace_file_on_interrupt)

# THIS PROGRAM IS FOR UNTRIMMED WHOLE ANALYSIS

bash_script = f"""#!/bin/bash

# Define the path of the potential trimmed file
TRIMMED_FILE="SRR_one/SRR_one_trimmed.fq.gz"

# Define the chromosomes to skip
SKIPPED_CHROMOSOMES="skip_chromosomes"

# Check if the trimmed file already exists
if [ ! -f "$TRIMMED_FILE" ]; then
    echo -e "\n\033[1;35mDownloading number sequence SRR_one from SRA...\033[0m "
    fastq-dump SRR_one

    if [ -d "SRR_one" ]; then
      rm -r "SRR_one"
    fi

    mkdir SRR_one
    mv SRR_one.fastq SRR_one

    echo -e "\n\033[1;35mRunning fastqc on SRR_one...\033[0m "
    fastqc SRR_one/SRR_one.fastq

    echo -e "\n\033[1;35mTrimming SRR_one...\033[0m "
    java -jar trim_path SE SRR_one/SRR_one.fastq SRR_one/SRR_one_trimmed.fq.gz ILLUMINACLIP:truseq3_path:2:30:10 SLIDINGWINDOW:4:20 MINLEN:35

    echo -e "\n\033[1;35mRunning fastqc on trimmed SRR_one...\033[0m "
    fastqc SRR_one/SRR_one_trimmed.fq.gz
else
    echo -e "\n\033[1;32mTrimmed file already exists. Skipping download, trimming, and quality check...\033[0m"
fi

echo -e "\n\033[1;35mMapping SRR_one reads using Bowtie2...\033[0m "
bowtie2 --very-fast-local -x bowtie_index_path SRR_one/SRR_one_trimmed.fq.gz -S SRR_one/SRR_one_mapped.sam

samtools view -S -b SRR_one/SRR_one_mapped.sam > SRR_one/SRR_one_mapped.bam

echo -e "\n\033[1;35mSorting using Samtools...\033[0m "
samtools sort SRR_one/SRR_one_mapped.bam > SRR_one/SRR_one_mapped.sorted.bam

# Check if the current chromosome should be skipped
if [[ ! $SKIPPED_CHROMOSOMES =~ $CURRENT_CHROMOSOME ]]; then
    echo -e "\n\033[1;35mSummarizing the base calls (mpileup)...\033[0m "
    bcftools mpileup -f bwa_chrom_path SRR_one/SRR_one_mapped.sorted.bam | bcftools call -mv -Ob -o SRR_one/SRR_one_mapped.raw.bcf

    echo -e "\n\033[1;35mFinalizing VCF...\033[0m "
    bcftools view SRR_one/SRR_one_mapped.raw.bcf | vcfutils.pl varFilter - > SRR_one/SRR_one_mapped.var.-final.vcf
fi

rm SRR_one/SRR_one.fastq SRR_one/SRR_one_mapped.sam SRR_one/SRR_one_mapped.bam SRR_one/SRR_one_mapped.sorted.bam SRR_one/SRR_one_mapped.raw.bcf SRR_one/SRR_one_fastqc.zip SRR_one/SRR_one_trimmed_fastqc.zip

"""

with open("untrimmed_bash_sra_v1.2.txt", "w") as f:
    f.write(bash_script)

# Define functions to replace text in files
def replace_text(file_path, old_text, new_text):
    with open(file_path, 'r+') as file:
        text = file.read().replace(old_text, new_text)
        file.seek(0)
        file.write(text)
        file.truncate()

def replace_in_untrimmed_bash_srr(old_text, new_text):
    replace_text('untrimmed_bash_sra_v1.2.txt', old_text, new_text)

# Get the current working directory
cwd = Path.cwd()
print('\n')

# Add the email to be notified when the process is done
user = input(f'{MAGENTA}1){RESET} Enter the email address to be notified once the analysis is complete: ')

# Add the job title
job_title = input(f'{MAGENTA}2){RESET} Enter a job name: ')

jar_file = '/usr/local/bin/Trimmomatic-0.39/trimmomatic-0.39.jar'
jar_path = os.path.expanduser(jar_file)

if os.path.exists(jar_path):
    trim_path = jar_path
    replace_in_untrimmed_bash_srr('trim_path', trim_path)
    print(f'{MAGENTA}3){RESET} {jar_path} is the absolute path.')
else:
    print(f'NOTE: {jar_path} does not match your absolute path.')
    print('You have a different path for trimmomatic-0.39.jar')
    trim_path = input(f'{MAGENTA}3){RESET} Copy and paste the absolute path to your trimmomatic-0.39.jar file: ')
    replace_in_untrimmed_bash_srr('trim_path', trim_path)

seq3_file = '/usr/local/bin/Trimmomatic-0.39/adapters/TruSeq3-SE.fa'
seq3_path = os.path.expanduser(seq3_file)

if os.path.exists(seq3_path):
    truseq3_path = seq3_path
    replace_in_untrimmed_bash_srr('truseq3_path', truseq3_path)
    print(f'{MAGENTA}4){RESET} {seq3_path} is the absolute path.')
else:
    print(f'NOTE: {seq3_path} does not match your absolute path.')
    print('You have a different path for TruSeq3-SE.fa')
    truseq3_path = input(f'{MAGENTA}4){RESET} Copy and paste the absolute path to your TruSeq3-SE.fa file: ')
    replace_in_untrimmed_bash_srr('truseq3_path', truseq3_path)

index_file = '/usr/local/bin/bowtie2-2.3.4.1-linux-x86_64/indexes/hg19'
index_path = os.path.expanduser(index_file)

if os.path.exists(index_path):
    bowtie_index_path = index_path
    replace_in_untrimmed_bash_srr('bowtie_index_path', bowtie_index_path)
    print(f'{MAGENTA}5){RESET} {index_path} is the absolute path.')
else:
    print(f'NOTE: {index_path} does not match your absolute path.')
    print('You have a different path for hg19')
    bowtie_index_path = input(f'{MAGENTA}5){RESET} Copy and paste the absolute path to your hg19 file: ')
    replace_in_untrimmed_bash_srr('bowtie_index_path', bowtie_index_path)

chrom_file = '/usr/local/bin/bwa.kit/resource-GRCh37/human_g1k_v37.fasta'
chrom_path = os.path.expanduser(chrom_file)

if os.path.exists(chrom_path):
    bwa_chrom_path = chrom_path
    replace_in_untrimmed_bash_srr('bwa_chrom_path', bwa_chrom_path)
    print(f'{MAGENTA}6){RESET} {chrom_path} is the absolute path.')
else:
    print(f'NOTE: {chrom_path} does not match your absolute path.')
    print('You have a different path for human_g1k_v37.fasta')
    bwa_chrom_path = input(f'{MAGENTA}6){RESET} Copy and paste the absolute path to your human_g1k_v37.fasta file: ')
    replace_in_untrimmed_bash_srr('bwa_chrom_path', bwa_chrom_path)

skip_chromosomes = input(f'{MAGENTA}7){RESET} Enter the chromosomes you would like to skip (separated by spaces): ')
replace_in_untrimmed_bash_srr('skip_chromosomes', skip_chromosomes)

# Save the paths and email to a file
with open('config.txt', 'w') as config:
    config.write(f'Email: {user}\n')
    config.write(f'Job Title: {job_title}\n')
    config.write(f'Trimmomatic path: {trim_path}\n')
    config.write(f'TruSeq3 path: {truseq3_path}\n')
    config.write(f'Bowtie2 index path: {bowtie_index_path}\n')
    config.write(f'BWA chromosome path: {bwa_chrom_path}\n')
    config.write(f'Skipped Chromosomes: {skip_chromosomes}\n')

# Make the bash script executable
os.chmod('untrimmed_bash_sra_v1.2.txt', 0o755)

# Rename the bash script
shutil.move('untrimmed_bash_sra_v1.2.txt', 'untrimmed_bash_sra_v1.2.sh')

print(f'{GREEN}\nConfiguration complete! You can now run the bash script "untrimmed_bash_sra_v1.2.sh".{RESET}')
