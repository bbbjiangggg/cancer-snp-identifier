#!/usr/bin/env python3

import os
import shutil
from pathlib import Path

# THIS PROGRAM IS FOR UNTRIMMED WHOLE ANALYSIS

# Define functions to replace text in files
def replace_text(file_path, old_text, new_text):
    with open(file_path, 'r+') as file:
        text = file.read().replace(old_text, new_text)
        file.seek(0)
        file.write(text)
        file.truncate()

def replace_in_untrimmed_bash_srr(old_text, new_text):
    replace_text('untrimmed_bash_srr.txt', old_text, new_text)


# Get the current working directory
cwd = Path.cwd()
print('\n')
print(f'\033[1;45m This is your current working directory: \033[0m{cwd}\n')

print('Use the command "realpath filename.txt" to get the complete path.\n')

# Add the email to be notified when the process is done
user = input('\033[1;45m 1) \033[0m Enter the email address to be notified once the analysis is complete: ')

# Add the job title
job_title = input('\033[1;45m 2) \033[0m Enter a job name: ')

# Add the path to where trimmomatic-0.39.jar is found
trim_path = input('\033[1;45m 3) \033[0m Copy and paste the complete path to your trimmomatic-0.39.jar file: ')
replace_in_untrimmed_bash_srr('trim_path', trim_path)

# Add the path to where TruSeq3 file is found
truseq3_path = input('\033[1;45m 4) \033[0m Copy and paste the complete path to your TruSeq3 file: ')
replace_in_untrimmed_bash_srr('truseq3_path', truseq3_path)

# Add the path to where bowtie files are found (must end in "bowtie")
bowtie2_path = input('\033[1;45m 5) \033[0m Copy and paste the complete path to your Bowtie files: ')
replace_in_untrimmed_bash_srr('bowtie2_path', bowtie2_path)

# Add the path to where reference chromosome is found
ref_chrom_path = input('\033[1;45m 6) \033[0m Copy and paste the complete path to your BWA reference chromosome: ')
replace_in_untrimmed_bash_srr('ref_chrom_path', ref_chrom_path)

# Make a copy of untrimmed bash srr
shutil.copy("untrimmed_bash_srr.txt", "copy_untrimmed_bash_srr.txt")


# This asks the user to type in the path to the accession list
accession = input("\033[1;45m 7) \033[0m Copy and paste the name of the accession list file (make sure it is in the same directory): ")
with open(accession, "r") as file:
    srr_list = [line.strip() for line in file if line.strip()]
    print(f"There are {len(srr_list)} unanalyzed sequences")
os.remove(f"copy_{accession}")


# This asks the user to type in the number of SRA sequences to be analyzed
num_sra_seqs = int(input("How many SRA sequences do you wish to analyze: "))

# This will store placement numbers into a list
ordinal = lambda n: f"{n}{['th', 'st', 'nd', 'rd', 'th', 'th', 'th', 'th', 'th', 'th'][n % 10 if n % 10 <= 3 and n % 100 not in (11, 12, 13) else 0]}"
placement = [ordinal(n) for n in range(1, num_sra_seqs + 1)]

# This will set different variables for different sra sequences
sra_list = srr_list[:num_sra_seqs]

# These commands will replace each SRA number on .txt file with each of the accession numbers entered by user
for index, sra in enumerate(sra_list):
    replace_in_untrimmed_bash_srr(f"accession_{index + 1}", sra)