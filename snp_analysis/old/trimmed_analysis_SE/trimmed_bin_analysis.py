#!/usr/bin/env python3

import os
import subprocess
from pathlib import Path
import signal
import importlib
import subprocess

packages = ['os', 'subprocess', 'pathlib', 'signal']

# Check if packages are installed, install them if necessary
for package in packages:
    try:
        importlib.import_module(package)
    except ImportError:
        print(f'{package} is not installed. Installing...')
        subprocess.run(['pip', 'install', package])

##############################################################################

# Define color codes
RED = '\033[1;31m'
GREEN = '\033[1;32m'
MAGENTA = '\033[1;35m'
BLUE = '\033[1;34m'
RESET = '\033[0m'

def replace_file_on_interrupt(sig, frame):
    # Remove the untrimmed_bash_sra_v1.2.txt file
    os.remove('trimmed_bash_sra_v1.2.txt')
    print('\nRemoving trimmed_bash_sra_v1.2.txt...')
    exit(1)

# Register the signal handler
signal.signal(signal.SIGINT, replace_file_on_interrupt)

################################################################################

# THIS PROGRAM IS FOR TRIMMED ANALYSIS

bash_script = f"""#!/bin/bash


bowtie2 --very-fast-local -x bowtie2_path srr_one/srr_one_trimmed.fq.gz -S srr_one/srr_one_mapped.sam

samtools view -S -b srr_one/srr_one_mapped.sam > srr_one/srr_one_mapped.bam

echo -e "\n\033[1;35mSorting using Samtools for srr_one...\033[0m "
samtools sort srr_one/srr_one_mapped.bam > srr_one/srr_one_mapped.sorted.bam

echo -e "\n\033[1;35mSummarizing the base calls (mpileup) for srr_one...\033[0m "
bcftools mpileup -f ref_chrom_path srr_one/srr_one_mapped.sorted.bam | bcftools call -mv -Ob -o srr_one/srr_one_mapped.raw.bcf

echo -e "\n\033[1;35mFinalizing VCF for srr_one...\033[0m "
bcftools view srr_one/srr_one_mapped.raw.bcf | vcfutils.pl varFilter - > srr_one/srr_one_mapped.var.-final.vcf

# Remove unnecessary files
rm srr_one/srr_one_mapped.sam srr_one/srr_one_mapped.bam srr_one/srr_one_mapped.sorted.bam srr_one/srr_one_mapped.raw.bcf

"""

with open("trimmed_bash_sra_v1.2.txt", "w") as f:
    f.write(bash_script)

##############################################################################

# THIS PROGRAM IS FOR TRIMMED ANALYSIS

# Define functions to replace text in files
def replace_text(file_path, old_text, new_text):
    with open(file_path, 'r+') as file:
        text = file.read().replace(old_text, new_text)
        file.seek(0)
        file.write(text)
        file.truncate()

def replace_in_trimmed_bash_sra(old_text, new_text):
    replace_text('trimmed_bash_sra_v1.2.txt', old_text, new_text)

# Get the current working directory
cwd = Path.cwd()
print('\n')
print(f'{MAGENTA} This is your current working directory: {RESET}{cwd}\n')

print("Use the command 'realpath filename.txt' to get the complete path.\n")

# Add the email to be notified when the process is done
user = input(f'{MAGENTA}1){RESET} Enter the email address to be notified once the analysis is complete: ')

# Add the job title
job_title = input(f'{MAGENTA}2){RESET} Enter a job name: ')


# Add the path to where bowtie files are found (must end in 'bowtie')
bowtie2_path = input(f'{MAGENTA}3){RESET} Copy and paste the complete path to your Bowtie files: ')
replace_in_trimmed_bash_sra('bowtie2_path', bowtie2_path)

# Add the path to where reference chromosome is found
ref_chrom_path = input(f'{MAGENTA}4){RESET} Copy and paste the complete path to your BWA reference chromosome: ')
replace_in_trimmed_bash_sra('ref_chrom_path', ref_chrom_path)


###########################################################################



# Define the directory to search in
search_dir = os.getcwd()

# Get a list of all directories in the search directory
src_dir = str(search_dir)


def find_unanalyzed_files(src_dir):
    files = [f for f in os.listdir(src_dir) if 'RR' in f and '.txt' not in f]

    vcfs = []
    for f in files:
        file_dir = os.path.join(src_dir, f)
        content = os.listdir(file_dir)
        for doc in content:
            if '_mapped.var.-final.vcf' in doc:
                name = doc.replace('_mapped.var.-final.vcf', '')
                vcfs.append(name)

    return sorted(list(set(files) - set(vcfs)))


def print_unanalyzed_files(unanalyzed_files, total_files):
    print(f'{MAGENTA}\nThese are the unanalyzed files in the current directory:{RESET}')
    for f in unanalyzed_files:
        print(f)
    print(f'{MAGENTA}There are{RESET} {len(unanalyzed_files)} {MAGENTA}unanalyzed files out of {total_files}.{RESET}')

# Find and print the sorted list of unanalyzed files
total_files = len(find_unanalyzed_files(search_dir))
new = find_unanalyzed_files(search_dir)
print_unanalyzed_files(new, total_files)

# This will store placement numbers into a list
ordinal = lambda n: "%d%s" % (n, "tsnrhtdd"[(n//10%10!=1)*(n%10<4)*n%10::4])
number = int(input('How many SRA sequences do you wish to analyze: '))
placement = [ordinal(n) for n in range(1, number + 1)]

# This will set different variables for different SRR sequences
srr_list = new[:number]

def replace_in_file(file_path, old_string, new_string):
    with open(file_path, 'r') as file:
        filedata = file.read()

    new_data = filedata.replace(old_string, new_string)

    with open(file_path, 'w') as file:
        file.write(new_data)


for index, srr in enumerate(srr_list):
    remaining = len(srr_list) - (index + 1)
    print(f"{MAGENTA}\nMapping reads using Bowtie2 for {srr}, {index+1} out of {len(srr_list)}, {remaining} left to analyze...{RESET}")
    replace_in_file('trimmed_bash_sra_v1.2.txt', 'number', placement[index])
    replace_in_file('trimmed_bash_sra_v1.2.txt', 'srr_one', srr)

    # Run the commands on the trimmed_bash_sra_v1.2.txt file
    os.system('cat trimmed_bash_sra_v1.2.txt | bash')

    # Replace the changed names back to the original
    replace_in_file('trimmed_bash_sra_v1.2.txt', placement[index], 'number')
    replace_in_file('trimmed_bash_sra_v1.2.txt', srr, 'srr_one')


    # Print how many accession numbers still need to be completed
    print(f"\n{srr} has been analyzed. {remaining} accession numbers still need to be completed.") 

# Reset the bowtie2_path and ref_chrom path
replace_in_file('trimmed_bash_sra_v1.2.txt', bowtie2_path, 'bowtie2_path')
replace_in_file('trimmed_bash_sra_v1.2.txt', ref_chrom_path, 'ref_chrom')





##############################################################################


#send an email to the user to let them ksrr_one the analysis is done
email_message = f"Your {job_title}_name analysis is complete. Please log in to check the results."
os.system(f'sendemail -f sudoroot1775@outlook.com -t {user} -u {job_title}_name Analysis Complete -m "{email_message}" -s smtp-mail.outlook.com:587 -o tls=yes -xu sudoroot1775@outlook.com -xp ydAEwVVu2s7uENC')
print(f'Sent email to {user}.')


def run_again():
    #runs trimmed_analysis_tools.py again
    print(f'{MAGENTA}\ntrimmed_analysis_tools.py is ready to run.{RESET}')

    while True:
        print(f'{MAGENTA}\nWould you like to run another analysis?{RESET}')
        choice = input('Enter yes or no to continue: ')

        if choice.lower() == 'yes':
            print(f'{MAGENTA}1){RESET} Email address used: ', user)
            print(f'{MAGENTA}2){RESET} Bowtie file path: ', bowtie2_path)
            print(f'{MAGENTA}3){RESET} BWA reference chromosome path: ', ref_chrom_path)
            
            # execute the analysis script
            os.system('python3 trimmed_analysis_tools_v1.2.py')
        
        elif choice.lower() == 'no':
            # exit the loop and end the program
            print(f'{MAGENTA}\nAnalysis terminated. Goodbye.{RESET}')
            break
        else:
            # handle invalid input by asking the user to enter either yes or no
            print(f'{MAGENTA}\nEnter either yes or no.{RESET}')

#Remove the trimmed_bash_sra_v1.2.txt file
os.system('rm trimmed_bash_sra_v1.2.txt')


run_again()