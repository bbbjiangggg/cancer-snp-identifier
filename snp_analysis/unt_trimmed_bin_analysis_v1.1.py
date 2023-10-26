#!/usr/bin/env python3

import os
import subprocess
from pathlib import Path
import signal

# Define color codes
RED = '\033[1;31m'
GREEN = '\033[1;32m'
MAGENTA = '\033[1;35m'
BLUE = '\033[1;34m'
RESET = '\033[0m'

#!/usr/bin/env python3

import os
import subprocess
from pathlib import Path
import signal

# Define color codes
RED = '\033[1;31m'
GREEN = '\033[1;32m'
MAGENTA = '\033[1;35m'
BLUE = '\033[1;34m'
RESET = '\033[0m'

def replace_file_on_interrupt(sig, frame):
    # Use absolute path for file operations
    bash_script_path = os.path.join(cwd, 'untrimmed_bash_sra_v1.2.txt')
    
    # Check if file exists before trying to remove it
    if os.path.exists(bash_script_path):
        os.remove(bash_script_path)
        print('\nRemoving untrimmed_bash_sra_v1.2.txt...')
    exit(1)

# Register the signal handler
signal.signal(signal.SIGINT, replace_file_on_interrupt)

# Define functions to replace text in files
def replace_text(file_path, old_text, new_text):
    # Check if the file exists before trying to open it
    if not os.path.exists(file_path):
        print(f"Error: The file {file_path} does not exist.")
        return
    try:
        with open(file_path, 'r+') as file:
            text = file.read().replace(old_text, new_text)
            file.seek(0)
            file.write(text)
            file.truncate()
    except Exception as e:
        print(f"An error occurred while trying to modify {file_path}: {str(e)}")


def replace_in_untrimmed_bash_srr(old_text, new_text):
    replace_text('untrimmed_bash_sra_v1.2.txt', old_text, new_text)

def get_absolute_path(prompt, default_path, variable_name):
    path = os.path.expanduser(default_path)
    if os.path.exists(path):
        print(f'{MAGENTA}{variable_name}){RESET} {path} is the absolute path.')
    else:
        print(f'NOTE: {path} does not match your absolute path.')
        print(f'You have a different path for {variable_name}')
        path = input(f'{MAGENTA}{variable_name}){RESET} Copy and paste the absolute path to your {variable_name} file: ')
    replace_in_untrimmed_bash_srr(variable_name, path)
    return path

def run_analysis(chromosome):
    # Update the paths in the bash script based on the chromosome
    bwa_chrom_path = f"/usr/local/bin/bwa/{chromosome}_bwa_ind/Homo_sapiens.GRCh38.dna.chromosome.{chromosome}.fa"
    bowtie_index_path = f"/usr/local/bin/bowtie/{chromosome}_bowtie_ind/bowtie"
    replace_in_untrimmed_bash_srr('bowtie_index_path', bowtie_index_path)
    replace_in_untrimmed_bash_srr('bwa_chrom_path', bwa_chrom_path)
    
    # Run the bash script
    subprocess.run(['bash', str(cwd) + '/untrimmed_bash_sra_v1.2.txt'])

    # Reset the paths in the bash script
    replace_in_untrimmed_bash_srr(bowtie_index_path, 'bowtie_index_path')
    replace_in_untrimmed_bash_srr(bwa_chrom_path, 'bwa_chrom_path')



# Define the paths
trim_path = get_absolute_path("Enter the absolute path to your trimmomatic-0.39.jar file:", '/usr/local/bin/Trimmomatic-0.39/trimmomatic-0.39.jar', 'trim_path')
truseq3_path = get_absolute_path("Enter the absolute path to your TruSeq3-SE.fa file:", '/usr/local/bin/Trimmomatic-0.39/adapters/TruSeq3-SE.fa', 'truseq3_path')

# Add the email to be notified when the process is done
user = input(f'{MAGENTA}1){RESET} Enter the email address to be notified once the analysis is complete: ')

# Add the job title
job_title = input(f'{MAGENTA}2){RESET} Enter a job name: ')

# This asks the user to type in the path to the accession list
accession = input(f'{MAGENTA}8){RESET} Copy and paste the name of the accession list file: ')

# Read the accession list file
with open(accession, 'r') as file:
    srr_list = [line.strip() for line in file if line.strip()]

# Create the bash script template
bash_script_template = """#!/bin/bash
# Define the path of the potential trimmed file
TRIMMED_FILE="SRR_one/SRR_one_trimmed.fq.gz"

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

echo -e "\n\033[1;35mSummarizing the base calls (mpileup)...\033[0m "
bcftools mpileup -f bwa_chrom_path SRR_one/SRR_one_mapped.sorted.bam | bcftools call -mv -Ob -o SRR_one/SRR_one_mapped.raw.bcf

echo -e "\n\033[1;35mFinalizing VCF...\033[0m "
bcftools view SRR_one/SRR_one_mapped.raw.bcf | vcfutils.pl varFilter - > SRR_one/SRR_one_mapped.var.-final.vcf

rm SRR_one/SRR_one.fastq SRR_one/SRR_one_mapped.sam SRR_one/SRR_one_mapped.bam SRR_one/SRR_one_mapped.sorted.bam SRR_one/SRR_one_mapped.raw.bcf SRR_one/SRR_one_fastqc.zip SRR_one/SRR_one_trimmed_fastqc.zip

"""

# Get the current working directory
cwd = Path.cwd()

# Create the bash script template file before proceeding
bash_script_path = os.path.join(cwd, 'untrimmed_bash_sra_v1.2.txt')
with open(bash_script_path, "w") as f:
    f.write(bash_script_template)

# Get the list of chromosomes to analyze
chromosomes = input(f"{MAGENTA}Enter the chromosome numbers (e.g., 1, 2, ... 22, X, Y) separated by commas to be analyzed: {RESET}").split(',')
chromosomes = [chrom.strip() for chrom in chromosomes]

# Create a new bash script for each chromosome, update paths, and run the analysis
for chromosome in chromosomes:
    # Update the paths in the bash script based on the chromosome
    bwa_chrom_path = f"/usr/local/bin/bwa/{chromosome}_bwa_ind/Homo_sapiens.GRCh38.dna.chromosome.{chromosome}.fa"
    bowtie_index_path = f"/usr/local/bin/bowtie/{chromosome}_bowtie_ind/bowtie"
    replace_in_untrimmed_bash_srr('bowtie_index_path', bowtie_index_path)
    replace_in_untrimmed_bash_srr('bwa_chrom_path', bwa_chrom_path)
    
    # Write the updated bash script to a file
    with open("untrimmed_bash_sra_v1.2.txt", "w") as f:
        f.write(bash_script_template)
    
    # Run the bash script
    subprocess.run(['bash', str(cwd) + '/untrimmed_bash_sra_v1.2.txt'])

    # Reset the paths in the bash script
    replace_in_untrimmed_bash_srr(bowtie_index_path, 'bowtie_index_path')
    replace_in_untrimmed_bash_srr(bwa_chrom_path, 'bwa_chrom_path')


# Send email notification
print('Sending email to ' + user + ' ....')
os.system('sendemail -f sudoroot1775@outlook.com -t ' + user + ' -u ' + job_title + '_Analysis Done -m "Ready to receive information for the next analysis." -s smtp-mail.outlook.com:587 -o tls=yes -xu sudoroot1775@outlook.com -xp ydAEwVVu2s7uENC')

# Clean up by removing the bash script file
os.remove('untrimmed_bash_sra_v1.2.txt')
