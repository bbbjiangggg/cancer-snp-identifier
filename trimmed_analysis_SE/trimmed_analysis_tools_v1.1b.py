#!/usr/bin/env python3
import os
import subprocess
from pathlib import Path

# THIS PROGRAM IS FOR TRIMMED FILES ONLY

# must have sendemail installed on terminal
# for Ubuntu use: $ sudo apt-get install libio-socket-ssl-perl libnet-ssleay-perl sendemail
# for Mac, use: brew install sendemail

# Define color codes
RED = '\033[1;31m'
GREEN = '\033[1;32m'
MAGENTA = '\033[1;35m'
BLUE = '\033[1;34m'
RESET = '\033[0m'

# Getting the current working directory
src_dir = os.getcwd()

# Printing current directory
print(f'{MAGENTA}\nCurrent working directory:{RESET} ' + src_dir)
print(f'{MAGENTA}\nUse the command "realpath file.txt" to get the complete path.{RESET}')

# Add the email to be notified when the process is done
user = input('Enter the email address to be notified once the analysis is complete: ')

# Add the job title
job = input('Enter a job name: ')

# Add the path to where bowtie files are found (must end in "bowtie/bowtie")
while True:
    bowtie = input('Copy and paste the complete path to your bowtie files: ')
    if os.path.isdir(bowtie) and os.path.isfile(os.path.join(bowtie, 'bowtie')):
        break
    else:
        print(f'{RED}Error: "{bowtie}" does not exist or is not a valid Bowtie directory. Please try again.{RESET}')

subprocess.run(['sed', '-i', f's|bowtie2_path/{bowtie}/g|{bowtie}|g', 'trimmed_bash_sra_v1.1.txt'])

# Add the path to where reference chromosome is found
while True:
    ref_chrom = input('Copy and paste the complete path to your reference chromosome: ')
    if os.path.isfile(ref_chrom):
        break
    else:
        print(f'{RED}Error: "{ref_chrom}" does not exist or is not a valid file. Please try again.{RESET}')

subprocess.run(['sed', '-i', f's|ref_chrom/{ref_chrom}/g', 'trimmed_bash_sra_v1.1.txt'])




# Get the path to the trimmed bash script
script_path = Path(__file__).parent / 'trimmed_bash_sra_v1.1.txt'

# Check if the script file exists
if not script_path.is_file():
    print(f"Error: Script file {script_path} not found")
else:
    # Replace the bowtie2_path variable in the script file
    subprocess.run(['sed', '-i', f's/bowtie2_path/{bowtie}/g', str(script_path)])

    # Replace the ref_chrom variable in the script file
    subprocess.run(['sed', '-i', f's/ref_chrom/{ref_chrom}/g', str(script_path)])

# Printing the sorted list of unanalyzed files
print(f'{MAGENTA}\nThese are the unanalyzed files in the current directory:{RESET}')

# This will store the list of unanalyzed files
files = os.listdir()
new = [i for i in files if 'RR' in i and '.txt' not in i]

vcfs = []
for j in new:
    srr_dirc = os.path.join(src_dir, j)
    content = os.listdir(srr_dirc)
    for doc in content:
        if '_mapped.var.-final.vcf' in doc:
            name = doc.replace('_mapped.var.-final.vcf', '')
            vcfs.append(name)

new = sorted(list(set(new) - set(vcfs)))
for i in new:
    print(i)

length = len(new)
print(f'*** There are {length} unanalyzed files ***')

# This will store placement numbers into a list
ordinal = lambda n: "%d%s" % (n, "tsnrhtdd"[(n//10%10!=1)*(n%10<4)*n%10::4])
number = int(input('How many SRA sequences do you wish to analyze: '))
placement = [ordinal(n) for n in range(1, number + 1)]

# This will set different variables for different SRR sequences
srr_list = new[:number]

# Replace the variables in the text file using subprocess
def replace_in_file(file_path, old_text, new_text):
    with open(file_path, 'r') as f:
        content = f.read()
    content = content.replace(old_text, new_text)
    with open(file_path, 'w') as f:
        f.write(content)

# These commands will replace each SRR number in the .txt file with 
# each of the accession numbers entered by the user
for index, srr in enumerate(srr_list):
    replace_in_file('trimmed_bash_sra_v1.1.txt', 'number', placement[index])
    replace_in_file('trimmed_bash_sra_v1.1.txt', 'now', srr)

    # Run the commands on the trimmed_bash_sra_v1.1.txt file
    command = 'cat trimmed_bash_sra_v1.1.txt | bash'
    subprocess.run(command, shell=True, check=True)

    # Replace the changed names back to the original
    replace_in_file('trimmed_bash_sra_v1.1.txt', placement[index], 'number')
    replace_in_file('trimmed_bash_sra_v1.1.txt', srr, 'now')

#reset the bowtie2_path and refchrome path
replace_in_file('trimmed_bash_sra_v1.1.txt', bowtie, 'bowtie2_path')
replace_in_file('trimmed_bash_sra_v1.1.txt', ref_chrom, 'ref_chrom')


#send an email to the user to let them know the analysis is done
email_message = f"Your {job}_name analysis is complete. Please log in to check the results."
sendemail_args = ['sendemail', '-f', 'sudoroot1775@outlook.com', '-t', user, '-u', f'{job}_name Analysis Complete', '-m', email_message, '-s', 'smtp-mail.outlook.com:587', '-o', 'tls=yes', '-xu', 'sudoroot1775@outlook.com', '-xp', 'ydAEwVVu2s7uENC']
subprocess.run(sendemail_args, check=True)

print(f'Sent email to {user}.')

print('\n')
print(f'{MAGENTA}\n trimmed_analysis_tools is ready to run.{RESET} \n')

def run_again():
    #runs trimmed_analysis_tools.py again
    print(f'{MAGENTA}\ntrimmed_analysis_tools.py is ready to run.{RESET} \n')

    while True:
        print(f'{MAGENTA}\nWould you like to run another analysis?\033[0;0;0m')
        choice = input('Enter yes or no to continue: ')

        if choice.lower() == 'yes':
            print(f'{MAGENTA}\n \033[0m Email address used:{RESET} ', user)
            print(f'{MAGENTA}\n \033[0m Bowtie file path:{RESET} ', bowtie)
            print(f'{MAGENTA}\n \033[0m BWA reference chromosome path:{RESET} ', ref_chrom)
            
            # execute the analysis script
            os.system('python3 trimmed_analysis_tools_v1.1.py')
        
        elif choice.lower() == 'no':
            # exit the loop and end the program
            print(f'{MAGENTA}\nAnalysis terminated. Goodbye.{RESET} \n')
            break
        else:
            # handle invalid input by asking the user to enter either yes or no
            print(f'{MAGENTA}\nEnter either yes or no.{RESET} \n')

run_again()