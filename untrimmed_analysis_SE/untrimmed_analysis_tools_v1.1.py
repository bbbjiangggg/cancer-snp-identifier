#!/usr/bin/env python3

import os
import shutil
import subprocess
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
    replace_text('untrimmed_bash_sra_v1.1.txt', old_text, new_text)


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
#trim_path = input('\033[1;45m 3) \033[0m Copy and paste the complete path to your trimmomatic-0.39.jar file: ')
#replace_in_untrimmed_bash_srr('trim_path', trim_path)




#trim_path = '~/Trimmomatic-0.39/trimmomatic-0.39.jar'
#replace_in_untrimmed_bash_srr('trim_path', trim_path)


jar_file = '~/Trimmomatic-0.39/trimmomatic-0.39.jar'
jar_path = os.path.expanduser(jar_file)


if os.path.exists(jar_path):
    trim_path = jar_path
    replace_in_untrimmed_bash_srr('trim_path', trim_path)
    print(f"\033[1;45m 3) \033[0m {jar_path} is the absolute path.")
else:
    print(f"NOTE: {jar_path} does not match your absolute path.")
    print("You have a different path for trimmomatic-0.39.jar")
    trim_path = input('\033[1;45m 3) \033[0m Copy and paste the absolute path to your trimmomatic-0.39.jar file: ')
    replace_in_untrimmed_bash_srr('trim_path', trim_path)



seq3_file = '~/Trimmomatic-0.39/adapters/TruSeq3-SE.fa'
seq3_path = os.path.expanduser(seq3_file)

if os.path.exists(seq3_path):
    truseq3_path = seq3_path
    replace_in_untrimmed_bash_srr('truseq3_path', truseq3_path)
    print(f"\033[1;45m 4) \033[0m {seq3_path} is the absolute path.")
else:
    print(f"NOTE: {seq3_path} does not match your absolute path.")
    print("You have a different path for TruSeq3-SE.fa")
    truseq3_path = input('\033[1;45m 4) \033[0m Copy and paste the absolute path to your TruSeq3-SE.fa file: ')
    replace_in_untrimmed_bash_srr('truseq3_path', truseq3_path)

# Add the path to where TruSeq3 file is found
#truseq3_path = input('\033[1;45m 4) \033[0m Copy and paste the complete path to your TruSeq3 file: ')
#replace_in_untrimmed_bash_srr('truseq3_path', truseq3_path)

# Add the path to where bowtie files are found (must end in "bowtie")
bowtie2_path = input('\033[1;45m 5) \033[0m Copy and paste the complete path to your Bowtie files: ')
replace_in_untrimmed_bash_srr('bowtie2_path', bowtie2_path)

# Add the path to where reference chromosome is found
ref_chrom_path = input('\033[1;45m 6) \033[0m Copy and paste the complete path to your BWA reference chromosome: ')
replace_in_untrimmed_bash_srr('ref_chrom_path', ref_chrom_path)

# Make a copy of untrimmed bash srr
shutil.copy("untrimmed_bash_sra_v1.1.txt", "copy_untrimmed_bash_sra_v1.1.txt")


# This asks the user to type in the path to the accession list
accession = input("\033[1;45m 7) \033[0m Copy and paste the name of the accession list file (make sure it is in the same directory): ")
with open(accession, "r") as file:
    srr_list = [line.strip() for line in file if line.strip()]
    print(f"There are {len(srr_list)} unanalyzed sequences")
os.remove("copy_untrimmed_bash_sra_v1.1.txt")


# This asks the user to type in the number of SRA sequences to be analyzed
num_sra_seqs = int(input("How many SRA sequences do you wish to analyze: "))

# This will store placement numbers into a list
ordinal = lambda n: f"{n}{['th', 'st', 'nd', 'rd', 'th', 'th', 'th', 'th', 'th', 'th'][n % 10 if n % 10 <= 3 and n % 100 not in (11, 12, 13) else 0]}"
placement = [ordinal(n) for n in range(1, num_sra_seqs + 1)]

# This will set different variables for different sra sequences
sra_list = srr_list[:num_sra_seqs]

# These commands will replace each SRA number on .txt file with each of the accession numbers entered by user
for index, sra in enumerate(sra_list):
    replace_in_untrimmed_bash_srr('number', placement[index])
    replace_in_untrimmed_bash_srr('SRR_one', sra)
    # Run the commands on the untrimmed_bash_sra_v1.1.txt file
    subprocess.run(['bash', str(cwd) + '/untrimmed_bash_sra_v1.1.txt'])

    # Replace the changed names back to original
    replace_in_untrimmed_bash_srr(placement[index], 'number')
    replace_in_untrimmed_bash_srr(sra, 'SRR_one')

# Reset the paths
replace_in_untrimmed_bash_srr(trim_path, 'trim_path')
replace_in_untrimmed_bash_srr(truseq3_path, 'truseq3_path')
replace_in_untrimmed_bash_srr(bowtie2_path, 'bowtie2_path')
replace_in_untrimmed_bash_srr(ref_chrom_path, 'ref_chrom_path')

# Remove the analyzed sequences from the accession list
with open(accession, "w") as file:
    for sra in srr_list[num_sra_seqs:]:
        file.write(sra + "\n")

# This will send an email to the user when the analysis is complete
#print(f'Sending email to {user}....')
#subprocess.run(['mail', '-s', f'{job_title} is complete', f'{user}'], input='Your analysis is complete', encoding='ascii')

#run the commands on the sendemail.txt file
print('Sending email to ' + user + ' ....')
os.system('sendemail -f sudoroot1775@outlook.com -t ' + user + ' -u ' + job_title + '_name_Analysis Done -m "Ready to receive information for the next analysis." -s smtp-mail.outlook.com:587 -o tls=yes -xu sudoroot1775@outlook.com -xp ydAEwVVu2s7uENC')


def run_analysis():
    """Runs untrimmed_analysis_tools.py and prompts user to continue or exit."""
    print('\033[1;45m untrimmed_analysis_tools is ready to run.\033[0;0;0m \n')
    while True:
        print('\033[1;45m Would you like to run another analysis? \033[0;0;0m')
        choice = input('Enter yes/no to continue: ')
        if choice.lower() == "yes":
            print('\033[1;45m 1) \033[0m Email address used: ', user)
            print('\033[1;45m 2) \033[0m trimmomatic-0.39.jar file path: ', trim_path)
            print('\033[1;45m 3) \033[0m TruSeq3 file path: ', truseq3_path)
            print('\033[1;45m 4) \033[0m Bowtie file path: ', bowtie2_path)
            print('\033[1;45m 5) \033[0m BWA reference chromosome path: ', ref_chrom_path)
            print('\033[1;45m 6) \033[0m Accession list file name: ', accession)
            os.system('python3 untrimmed_analysis_tools_v1.1.py')

        elif choice.lower() == "no":
            print('\033[1;45m Analysis terminated. Goodbye. \033[0;0;0m')
            break
        else:
            print('\033[1;45m Enter either yes or no.\033[0;0;0m ')

run_analysis()

