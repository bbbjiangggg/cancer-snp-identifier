
import os
import subprocess
import signal
from pathlib import Path
import snp_analysis

# Define color codes
RED = '\033[1;31m'
GREEN = '\033[1;32m'
MAGENTA = '\033[1;35m'
BLUE = '\033[1;34m'
RESET = '\033[0m'

def replace_file_on_interrupt(sig, frame):
    print('\nProcess interrupted by user...')
    exit(1)

# Register the signal handler
signal.signal(signal.SIGINT, replace_file_on_interrupt)

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

# Initialize paths
trim_path = "/usr/local/bin/Trimmomatic-0.39/trimmomatic-0.39.jar"
truseq3_path = "/usr/local/bin/Trimmomatic-0.39/adapters/TruSeq3-SE.fa"

# Ask the user if they want to analyze the whole genome or a specific chromosome
analysis_scope = input(f"{MAGENTA}Would you like to analyze the whole genome or a specific chromosome?\n1) Whole Genome\n2) Specific Chromosome\nEnter the number corresponding to your choice: {RESET}")
bwa_chrom_path = "/usr/local/bin/bwa/hg38/GRCh38_reference.fa" if analysis_scope == "1" else "/usr/local/bin/bwa/chrom_bwa_ind/Homo_sapiens.GRCh38.dna.chromosome.chrom.fa"
bowtie_index_path = "/usr/local/bin/bowtie/hg38/bowtie" if analysis_scope == "1" else "/usr/local/bin/bowtie/chrom_bowtie_ind/bowtie"

# This asks the user to type in the path to the accession list
accession = input(f'{MAGENTA}8){RESET} Copy and paste the name of the accession list file: ')

# Read the accession list file
with open(accession, 'r') as file:
    srr_list = [line.strip() for line in file if line.strip()]

# Create a new list of accession numbers to be analyzed
to_analyze = []
for sra in srr_list:
    # Check if a directory with the same SRA/ERR number exists and has a file ending with mapped.var.-final.vcf
    sra_dir = f'{cwd}/{sra}'
    vcf_file = f'{sra_dir}/{sra}_mapped.var.-final.vcf'
    if os.path.exists(sra_dir) and os.path.isfile(vcf_file):
        print(f'{sra} already has a mapped.var.-final.vcf file in the directory. Skipping analysis...')
    else:
        to_analyze.append(sra)

# This asks the user to type in the number of SRA sequences to be analyzed
num_sra_seqs = int(input(f'{MAGENTA}9){RESET}How many SRA sequences do you wish to analyze (out of {len(to_analyze)} remaining)? '))

# Set different variables for different sra sequences
sra_list = to_analyze[:num_sra_seqs]

# Run the SNP analysis
for sra in sra_list:
    snp_analysis.run_analysis(
        sra=sra,
        trim_path=trim_path,
        truseq3_path=truseq3_path,
        bowtie_index_path=bowtie_index_path,
        bwa_chrom_path=bwa_chrom_path,
        output_dir=cwd
    )

# Send email notification
print('Sending email to ' + user + ' ....')
os.system('sendemail -f sudoroot1775@outlook.com -t ' + user + ' -u ' + job_title + '_Analysis Done -m "Ready to receive information for the next analysis." -s smtp-mail.outlook.com:587 -o tls=yes -xu sudoroot1775@outlook.com -xp ydAEwVVu2s7uENC')

def run_analysis():
    # Runs untrimmed_analysis_tools.py and prompts user to continue or exit.
    print(f'{MAGENTA}untrimmed_analysis_tools is ready to run.{RESET} \n')
    while True:
        print(f'{MAGENTA}Would you like to run another analysis? {RESET}')
        choice = input('Enter yes/no to continue: ')
        if choice.lower() == 'yes':
            print(f'{MAGENTA} 1) {RESET} Email address used: ', user)
            print(f'{MAGENTA} 2) {RESET} trimmomatic-0.39.jar file path: ', trim_path)
            print(f'{MAGENTA} 3) {RESET} TruSeq3 file path: ', truseq3_path)
            print(f'{MAGENTA} 4) {RESET} Bowtie file path: ', bowtie_index_path)
            print(f'{MAGENTA} 5) {RESET} BWA reference chromosome path: ', bwa_chrom_path)
            print(f'{MAGENTA} 6) {RESET} Accession list file name: ', accession)
            os.system('python3 /usr/local/bin/untrimmed_bin_analysis.py')

        elif choice.lower() == 'no':
            print(f'{MAGENTA} Analysis terminated. Goodbye. {RESET}')
            break
        else:
            print(f'{MAGENTA} Enter either yes or no.{RESET} ')

run_analysis()
