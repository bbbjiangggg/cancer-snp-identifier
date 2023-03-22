#!/usr/bin/env python3

import os
import shutil
import subprocess
from pathlib import Path
import argparse

# THIS PROGRAM IS FOR UNTRIMMED WHOLE ANALYSIS

# Define color codes
RED = '\033[1;31m'
GREEN = '\033[1;32m'
MAGENTA = '\033[1;35m'
BLUE = '\033[1;34m'
RESET = '\033[0m'

# Define functions to replace text in files
def replace_text(file_path, old_text, new_text):
    with open(file_path, 'r+') as file:
        text = file.read().replace(old_text, new_text)
        file.seek(0)
        file.write(text)
        file.truncate()

def replace_in_untrimmed_bash_srr(old_text, new_text):
    replace_text('untrimmed_bash_sra_v1.2.txt', old_text, new_text)

# Define the function to find missing VCF directories
def find_missing_vcf_directories(path="."):
    missing_vcf_dirs = []

    for dir_name in os.listdir(path):
        if dir_name.startswith("SRR") or dir_name.startswith("ERR"):
            dir_path = os.path.join(path, dir_name)
            if os.path.isdir(dir_path):
                vcf_files = [file for file in os.listdir(dir_path) if file.endswith("_mapped.var.-final.vcf")]
                if not vcf_files:
                    missing_vcf_dirs.append(dir_path)

    return missing_vcf_dirs


# Get the current working directory
cwd = Path.cwd()
print('\n')
print(f'This is your current working directory: {cwd}')

print('Use the command "realpath filename.txt" to get the complete path.\n')

# Create argument parser
parser = argparse.ArgumentParser(description='Untrimmed analysis tools')
parser.add_argument('email', help='email address to be notified once the analysis is complete')
parser.add_argument('job_title', help='job name')
parser.add_argument('--trimmomatic', help='absolute path to trimmomatic-0.39.jar file')
parser.add_argument('--truseq3', help='absolute path to TruSeq3-SE.fa file')
parser.add_argument('--bowtie', help='absolute path to Bowtie files')
parser.add_argument('--ref', help='absolute path to BWA reference chromosome')
parser.add_argument('--accession', help='name of the accession list file')
parser.add_argument('--num-sra-seqs', type=int, default=1, help='number of SRA sequences to analyze')

# Parse arguments
args = parser.parse_args()

# Add the email to be notified when the process is done
email_me = input(f'{MAGENTA}1){RESET} Enter your email address to be notified when the analysis is complete: ')
replace_in_untrimmed_bash_srr('email_me', email_me)
user = args.email

# Add the job title
job_title = args.job_title

jar_file = '~/Trimmomatic-0.39/trimmomatic-0.39.jar'
if args.trimmomatic:
    jar_file = args.trimmomatic
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

    seq3_file = '~/Trimmomatic-0.39/adapters/TruSeq3-SE.fa'
    seq3_path = Path(seq3_file).expanduser()

    if seq3_path.exists():
        truseq3_path = seq3_path
        replace_in_untrimmed_bash_srr('truseq3_path', truseq3_path)
        print(f'{MAGENTA}4){RESET} {seq3_path} is the absolute path.')
    else:
        print(f'NOTE: {seq3_path} does not match your absolute path.')
        print('You have a different path for TruSeq3-SE.fa')
        truseq3_path = input(f'{MAGENTA}4){RESET} Copy and paste the absolute path to your TruSeq3-SE.fa file: ')
        replace_in_untrimmed_bash_srr('truseq3_path', truseq3_path)

# Add the path to where bowtie files are found (must end in "bowtie")
bowtie2_path = args.bowtie_path
replace_in_untrimmed_bash_srr('bowtie2_path', bowtie2_path)

# Add the path to where reference chromosome is found
ref_chrom_path = args.ref_chrom_path
replace_in_untrimmed_bash_srr('ref_chrom_path', ref_chrom_path)

# Make a copy of untrimmed bash srr
shutil.copy("untrimmed_bash_sra_v1.2.txt", "copy_untrimmed_bash_sra_v1.2.txt")

# Get the list of directories missing the required VCF file
missing_vcf_dirs = find_missing_vcf_directories(args.path_to_check)

# This asks the user to type in the path to the accession list
accession = args.accession_file
with open(accession, "r") as file:
    srr_list = [line.strip() for line in file if line.strip()]
    print(f'There are {len(srr_list)} unanalyzed sequences')
os.remove('copy_untrimmed_bash_sra_v1.2.txt')

# This asks the user to type in the number of SRA sequences to be analyzed
num_sra_seqs = args.num_sra_seqs

# This will store placement numbers into a list
ordinal = lambda n: f"{n}{['th', 'st', 'nd', 'rd', 'th', 'th', 'th', 'th', 'th', 'th'][n % 10 if n % 10 <= 3 and n % 100 not in (11, 12, 13) else 0]}"
placement = [ordinal(n) for n in range(1, num_sra_seqs + 1)]

# This will set different variables for different sra sequences
sra_list = srr_list[:num_sra_seqs]

# These commands will replace each SRA number on .txt file with each of the accession numbers entered by user
for index, sra in enumerate(sra_list):
    replace_in_untrimmed_bash_srr('number', placement[index])
    replace_in_untrimmed_bash_srr('SRR_one', sra)
    # Run the commands on the untrimmed_bash_sra_v1.2.txt file
    subprocess.run(['bash', str(cwd) + '/untrimmed_bash_sra_v1.2.txt'])

    # Replace the changed names back to original
    replace_in_untrimmed_bash_srr(placement[index], 'number')
    replace_in_untrimmed_bash_srr(sra, 'SRR_one')

# Reset the paths
replace_in_untrimmed_bash_srr(trim_path, 'trim_path')
replace_in_untrimmed_bash_srr(truseq3_path, 'truseq3_path')
replace_in_untrimmed_bash_srr(bowtie2_path, 'bowtie2_path')
replace_in_untrimmed_bash_srr(ref_chrom_path, 'ref_chrom_path')
replace_in_untrimmed_bash_srr(placement[index], 'number')
replace_in_untrimmed_bash_srr(sra, 'SRR_one')

# Remove the analyzed sequences from the accession list
'''with open(accession, "w") as file:
    for sra in srr_list[num_sra_seqs:]:
        file.write(sra + "\n")'''

# This will send an email to the user when the analysis is complete
#print(f'Sending email to {user}....')
#subprocess.run(['mail', '-s', f'{job_title} is complete', f'{user}'], input='Your analysis is complete', encoding='ascii')

#run the commands on the sendemail.txt file
print('Sending email to ' + user + ' ....')
os.system('sendemail -f sudoroot1775@outlook.com -t ' + user + ' -u ' + job_title + '_name_Analysis Done -m "Ready to receive information for the next analysis." -s smtp-mail.outlook.com:587 -o tls=yes -xu sudoroot1775@outlook.com -xp ydAEwVVu2s7uENC')


def run_analysis():
    #Runs untrimmed_analysis_tools.py and prompts user to continue or exit.
    print(f'{MAGENTA}\n untrimmed_analysis_tools is ready to run.{RESET}')
    while True:
        print(f'{MAGENTA}\n Would you like to run another analysis? {RESET}')
        choice = input('Enter yes/no to continue: ')
        if choice.lower() == "yes":
            print(f'{MAGENTA} 1) {RESET} Email address used: ', user)
            print(f'{MAGENTA} 2) {RESET} trimmomatic-0.39.jar file path: ', trim_path)
            print(f'{MAGENTA} 3) {RESET} TruSeq3 file path: ', truseq3_path)
            print(f'{MAGENTA} 4) {RESET} Bowtie file path: ', bowtie2_path)
            print(f'{MAGENTA} 5) {RESET} BWA reference chromosome path: ', ref_chrom_path)
            print(f'{MAGENTA} 6) {RESET} Accession list file name: ', accession)
            os.system('python3 untrimmed_analysis_tools_v1.1.py')

        elif choice.lower() == "no":
            print(f'{MAGENTA}\n Analysis terminated. Goodbye. {RESET}')
            break
        else:
            print(f'{MAGENTA}\n Enter either yes or no.{RESET} ')

run_analysis()

